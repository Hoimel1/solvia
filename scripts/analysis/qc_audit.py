#!/usr/bin/env python3
"""Aggregate PMF QC results and flag runs that need resimulation.

This script scans simulation run directories (default: ``simulations``),
parses the existing ``qc_report.yaml`` and ``pmf_analysis_results.yaml``
artifacts, and generates human-friendly summaries (CSV + Markdown).

- Low effective sample size (ESS) windows, insufficient overlap, replicate
  conflicts, and tilted bulk profiles are surfaced explicitly.
- Runs crossing the hard thresholds are marked as ``resimulate``; optional
  ``--mark-resim`` creates a ``RESIMULATE.flag`` marker next to the PMF tag.

The script is importable for testing; invoke ``main()`` only when executed
directly from the CLI.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import warnings
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import yaml


@dataclass(frozen=True)
class QCThresholds:
    """Numeric thresholds for QC decisions."""

    ess_fail: int = 150
    ess_warn: int = 200
    overlap_fail: float = 0.15
    tilt_span: float = 12.0
    bulk_start_nm: float = 2.2
    bulk_min_points: int = 20


@dataclass
class QCIssue:
    """Single QC issue raised for a run/tag."""

    type: str
    severity: str
    message: str
    data: Dict[str, object] = field(default_factory=dict)


@dataclass
class RunQC:
    """Aggregated QC metrics for a single run + PMF tag."""

    run_dir: Path
    tag_dir: Path
    run_name: str
    tag: str
    replicate: Optional[int]
    peptide_id: str
    min_ess: Optional[float] = None
    required_ess: Optional[float] = None
    min_overlap: Optional[float] = None
    ess_fail_centers: List[float] = field(default_factory=list)
    overlap_fail_pairs: List[Tuple[float, float]] = field(default_factory=list)
    replicate_conflicts: List[Dict[str, object]] = field(default_factory=list)
    tilt_span: Optional[float] = None
    tilt_slope: Optional[float] = None
    issues: List[QCIssue] = field(default_factory=list)

    @property
    def status(self) -> str:
        has_errors = any(issue.severity == "error" for issue in self.issues)
        if has_errors:
            return "resimulate"
        if self.issues:
            return "review"
        return "ok"

    @property
    def reasons(self) -> List[str]:
        return [issue.type for issue in self.issues]


def _safe_load_yaml(path: Path) -> Optional[dict]:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _extract_run_metadata(run_dir: Path) -> Tuple[str, Optional[int]]:
    """Return (peptide_id, replicate) parsed from the run directory name."""

    stem = run_dir.name
    # Expected format: solvia_<id>_run_<replicate>
    parts = stem.split("_run_")
    if len(parts) == 2:
        peptide = parts[0]
        try:
            replicate = int(parts[1])
        except ValueError:
            replicate = None
        return peptide.upper(), replicate
    return stem.upper(), None


def _analyze_qc_report(run: RunQC, report: dict, thresholds: QCThresholds) -> None:
    ess_entries = [
        entry
        for entry in report.get("ess_check", [])
        if not entry.get("superseded", False)
    ]
    if ess_entries:
        run.min_ess = float(min(entry.get("ess", math.inf) for entry in ess_entries))
        # min_required reported per entry; take minimum for context
        required_vals = [entry.get("min_required") for entry in ess_entries]
        required_vals = [float(v) for v in required_vals if v is not None]
        run.required_ess = min(required_vals) if required_vals else None
        for entry in ess_entries:
            if not bool(entry.get("passed", False)):
                run.ess_fail_centers.append(float(entry.get("center")))
        if run.ess_fail_centers:
            run.issues.append(
                QCIssue(
                    type="low_ess",
                    severity="error",
                    message="ESS below threshold",
                    data={"centers": run.ess_fail_centers},
                )
            )
        elif run.min_ess is not None and run.min_ess < thresholds.ess_warn:
            run.issues.append(
                QCIssue(
                    type="ess_margin",
                    severity="warn",
                    message=f"Low ESS headroom (min {run.min_ess:.1f})",
                    data={},
                )
            )

    overlap_entries = report.get("overlap_check", [])
    if overlap_entries:
        # Some entries may be duplicated for superseded replicates; retain minimum
        mins = []
        for entry in overlap_entries:
            val = entry.get("overlap")
            if val is None:
                continue
            try:
                mins.append(float(val))
            except (TypeError, ValueError):
                continue
            if not bool(entry.get("passed", False)):
                wins = entry.get("windows") or []
                if len(wins) == 2:
                    run.overlap_fail_pairs.append((float(wins[0]), float(wins[1])))
        run.min_overlap = min(mins) if mins else None
        if run.overlap_fail_pairs:
            run.issues.append(
                QCIssue(
                    type="low_overlap",
                    severity="error",
                    message="Adjacent window overlap below gate",
                    data={"pairs": run.overlap_fail_pairs},
                )
            )
        elif run.min_overlap is not None and run.min_overlap < thresholds.overlap_fail:
            run.issues.append(
                QCIssue(
                    type="overlap_margin",
                    severity="warn",
                    message=f"Overlap close to limit (min {run.min_overlap:.3f})",
                    data={},
                )
            )


def _detect_replicate_conflicts(run: RunQC, metadata: dict) -> None:
    selection = metadata.get("qc_selection")
    if not isinstance(selection, dict):
        return
    conflicts: List[Dict[str, object]] = []
    for center_raw, info in selection.items():
        reps = info.get("replicates") or []
        drops = [rep for rep in reps if rep.get("drop")]
        reasons = []
        for rep in reps:
            if rep.get("reason"):
                reasons.extend(rep["reason"])
        if drops or reasons:
            try:
                center = float(center_raw)
            except (TypeError, ValueError):
                center = center_raw
            conflicts.append({
                "center": center,
                "dropped": len(drops),
                "reasons": reasons,
            })
    if conflicts:
        run.replicate_conflicts = conflicts
        run.issues.append(
            QCIssue(
                type="replicate_conflict",
                severity="error",
                message="Replicates filtered due to QC",
                data={"conflicts": conflicts},
            )
        )


def _compute_tilt_metrics(run: RunQC, analysis: dict, thresholds: QCThresholds) -> None:
    pmf_data = analysis.get("pmf_data") if isinstance(analysis, dict) else None
    if not isinstance(pmf_data, dict):
        return
    z_vals = pmf_data.get("z")
    pmf_vals = pmf_data.get("pmf")
    if not (isinstance(z_vals, Sequence) and isinstance(pmf_vals, Sequence)):
        return
    z = np.asarray(z_vals, dtype=float)
    pmf = np.asarray(pmf_vals, dtype=float)
    mask = (~np.isnan(pmf)) & (z >= thresholds.bulk_start_nm)
    if int(np.sum(mask)) < thresholds.bulk_min_points:
        return
    bulk = pmf[mask]
    run.tilt_span = float(np.nanmax(bulk) - np.nanmin(bulk))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", np.RankWarning)
        slope = np.polyfit(z[mask], bulk, 1)[0]
    run.tilt_slope = float(slope)
    if run.tilt_span is not None and run.tilt_span >= thresholds.tilt_span:
        run.issues.append(
            QCIssue(
                type="tilted_pmf",
                severity="error",
                message=f"Bulk span {run.tilt_span:.1f} kJ/mol exceeds limit",
                data={"span": run.tilt_span, "slope": run.tilt_slope},
            )
        )


def analyze_simulation_runs(
    sim_root: Path,
    thresholds: QCThresholds,
) -> List[RunQC]:
    """Inspect all run directories and collect QC summaries."""

    results: List[RunQC] = []
    if not sim_root.exists():
        return results

    for run_dir in sorted(sim_root.glob("solvia_*_run_*")):
        pmf_root = run_dir / "pmf"
        if not pmf_root.exists():
            continue
        peptide_id, replicate = _extract_run_metadata(run_dir)
        for tag_dir in sorted(pmf_root.iterdir()):
            if not tag_dir.is_dir():
                continue
            qc_file = tag_dir / "qc_report.yaml"
            metadata_file = tag_dir / "pmf_metadata.yaml"
            analysis_file = tag_dir / "pmf_analysis_results.yaml"
            if not qc_file.exists():
                continue
            run = RunQC(
                run_dir=run_dir,
                tag_dir=tag_dir,
                run_name=run_dir.name,
                tag=tag_dir.name,
                replicate=replicate,
                peptide_id=peptide_id,
            )
            report = _safe_load_yaml(qc_file) or {}
            _analyze_qc_report(run, report, thresholds)
            metadata = _safe_load_yaml(metadata_file) or {}
            _detect_replicate_conflicts(run, metadata)
            analysis = _safe_load_yaml(analysis_file)
            if analysis:
                _compute_tilt_metrics(run, analysis, thresholds)
            results.append(run)
    return results


def _ensure_output_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_reports(
    runs: Sequence[RunQC],
    output_dir: Path,
    mark_resim: bool = False,
) -> Path:
    """Write CSV, Markdown, and YAML manifests summarising QC findings."""

    _ensure_output_dir(output_dir)
    csv_path = output_dir / "qc_report.csv"
    md_path = output_dir / "qc_report.md"
    yaml_path = output_dir / "qc_resimulate.yaml"

    fieldnames = [
        "run",
        "tag",
        "status",
        "reasons",
        "min_ess",
        "required_ess",
        "ess_fail_centers",
        "min_overlap",
        "overlap_fail_pairs",
        "replicate_conflicts",
        "tilt_span",
        "tilt_slope",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for run in runs:
            writer.writerow(
                {
                    "run": run.run_name,
                    "tag": run.tag,
                    "status": run.status,
                    "reasons": ",".join(run.reasons),
                    "min_ess": f"{run.min_ess:.1f}" if run.min_ess is not None else "",
                    "required_ess": f"{run.required_ess:.1f}" if run.required_ess else "",
                    "ess_fail_centers": ";".join(
                        f"{c:.2f}" for c in run.ess_fail_centers
                    ),
                    "min_overlap": f"{run.min_overlap:.3f}" if run.min_overlap is not None else "",
                    "overlap_fail_pairs": ";".join(
                        f"({a:.2f},{b:.2f})" for a, b in run.overlap_fail_pairs
                    ),
                    "replicate_conflicts": ";".join(
                        f"z={c['center']}\u2192{c['dropped']}" for c in run.replicate_conflicts
                    ),
                    "tilt_span": f"{run.tilt_span:.1f}" if run.tilt_span is not None else "",
                    "tilt_slope": f"{run.tilt_slope:.2f}" if run.tilt_slope is not None else "",
                }
            )

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M %Z")
    status_counts: Dict[str, int] = {"resimulate": 0, "review": 0, "ok": 0}
    for run in runs:
        status_counts[run.status] = status_counts.get(run.status, 0) + 1

    lines: List[str] = []
    lines.append(f"# QC Audit Report\n")
    lines.append(f"Generated: {timestamp}\n")
    lines.append("## Summary\n")
    lines.append("Status | Count\n")
    lines.append("--- | ---\n")
    for status in ("resimulate", "review", "ok"):
        lines.append(f"{status} | {status_counts.get(status, 0)}\n")
    lines.append("\n")

    def _section(title: str, candidates: Iterable[RunQC]) -> None:
        items = list(candidates)
        if not items:
            return
        lines.append(f"## {title}\n")
        for run in items:
            reason_txt = ", ".join(run.reasons) or "(none)"
            detail_parts = []
            if run.min_ess is not None:
                detail_parts.append(f"min_ESS={run.min_ess:.1f}")
            if run.min_overlap is not None:
                detail_parts.append(f"min_overlap={run.min_overlap:.3f}")
            if run.tilt_span is not None:
                detail_parts.append(f"tilt_span={run.tilt_span:.1f}")
            lines.append(
                f"- `{run.run_name}/{run.tag}` → {run.status} ({reason_txt}; "
                + ", ".join(detail_parts)
                + ")\n"
            )
        lines.append("\n")

    _section("Needs Resimulation", (r for r in runs if r.status == "resimulate"))
    _section("Review Manually", (r for r in runs if r.status == "review"))

    md_path.write_text("".join(lines), encoding="utf-8")

    resim_entries = [
        {
            "run": run.run_name,
            "tag": run.tag,
            "reasons": run.reasons,
        }
        for run in runs
        if run.status == "resimulate"
    ]
    with yaml_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump({"resimulate": resim_entries}, handle, sort_keys=False)

    if mark_resim:
        for run in runs:
            if run.status != "resimulate":
                continue
            flag_path = run.tag_dir / "RESIMULATE.flag"
            message_lines = [
                "qc_audit: resimulation required",
                f"reasons: {', '.join(run.reasons)}",
                f"generated: {timestamp}",
            ]
            flag_path.write_text("\n".join(message_lines) + "\n", encoding="utf-8")

    return md_path


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate PMF QC diagnostics")
    parser.add_argument(
        "sim_root",
        nargs="?",
        type=Path,
        default=Path("simulations"),
        help="Root directory containing solvia_* run folders",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis/qc"),
        help="Where to write qc_report.csv/md",
    )
    parser.add_argument(
        "--tilt-span",
        type=float,
        default=QCThresholds.tilt_span,
        help="Tilt span threshold in kJ/mol for flagging runs",
    )
    parser.add_argument(
        "--ess-fail",
        type=int,
        default=QCThresholds.ess_fail,
        help="ESS threshold regarded as a hard failure",
    )
    parser.add_argument(
        "--ess-warn",
        type=int,
        default=QCThresholds.ess_warn,
        help="ESS threshold for warnings",
    )
    parser.add_argument(
        "--overlap-min",
        type=float,
        default=QCThresholds.overlap_fail,
        help="Minimum acceptable adjacent overlap",
    )
    parser.add_argument(
        "--bulk-start",
        type=float,
        default=QCThresholds.bulk_start_nm,
        help="Lower bound (nm) for bulk tilt evaluation",
    )
    parser.add_argument(
        "--bulk-min-points",
        type=int,
        default=QCThresholds.bulk_min_points,
        help="Minimum histogram points required for tilt analysis",
    )
    parser.add_argument(
        "--mark-resim",
        action="store_true",
        help="Create RESIMULATE.flag files for failing runs",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    thresholds = QCThresholds(
        ess_fail=args.ess_fail,
        ess_warn=args.ess_warn,
        overlap_fail=args.overlap_min,
        tilt_span=args.tilt_span,
        bulk_start_nm=args.bulk_start,
        bulk_min_points=args.bulk_min_points,
    )
    runs = analyze_simulation_runs(args.sim_root, thresholds)
    if not runs:
        print(f"No runs found under {args.sim_root}")
        return 0
    write_reports(runs, args.output_dir, mark_resim=args.mark_resim)
    print(f"QC report written to {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
