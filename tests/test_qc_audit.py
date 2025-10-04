import importlib.util
import sys
from pathlib import Path


def load_module():
    module_path = Path("scripts/analysis/qc_audit.py")
    spec = importlib.util.spec_from_file_location("qc_audit", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def make_run(tmp_path: Path) -> Path:
    run_dir = tmp_path / "solvia_999_run_1" / "pmf" / "pmf_midplane"
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def write_minimal_files(run_dir: Path) -> None:
    qc_report = {
        "ess_check": [
            {
                "center": 0.0,
                "ess": 120,
                "min_required": 150,
                "passed": False,
                "superseded": False,
            }
        ],
        "overlap_check": [
            {
                "windows": [0.0, 0.1],
                "overlap": 0.25,
                "passed": True,
            }
        ],
    }
    metadata = {
        "qc_selection": {
            0.0: {
                "replicates": [
                    {
                        "center": 0.0,
                        "drop": True,
                        "reason": ["low_ess"],
                    }
                ]
            }
        }
    }
    analysis = {
        "pmf_data": {
            "z": [0.0, 2.2, 2.3, 2.4, 2.5],
            "pmf": [0.0, -8.0, -6.0, -4.0, 6.0],
        }
    }
    import yaml

    with (run_dir / "qc_report.yaml").open("w") as handle:
        yaml.safe_dump(qc_report, handle)
    with (run_dir / "pmf_metadata.yaml").open("w") as handle:
        yaml.safe_dump(metadata, handle)
    with (run_dir / "pmf_analysis_results.yaml").open("w") as handle:
        yaml.safe_dump(analysis, handle)


def test_qc_audit_flags_low_ess_and_tilt(tmp_path):
    mod = load_module()
    run_dir = make_run(tmp_path)
    write_minimal_files(run_dir)

    thresholds = mod.QCThresholds(
        ess_fail=150,
        ess_warn=180,
        overlap_fail=0.2,
        tilt_span=5.0,
        bulk_start_nm=2.0,
        bulk_min_points=3,
    )
    results = mod.analyze_simulation_runs(tmp_path, thresholds)
    assert len(results) == 1
    summary = results[0]
    assert summary.status == "resimulate"
    assert {"low_ess", "replicate_conflict", "tilted_pmf"}.issubset(summary.reasons)
    assert summary.min_ess == 120.0
    assert summary.replicate_conflicts


def test_write_reports_creates_outputs(tmp_path):
    mod = load_module()
    run_dir = make_run(tmp_path)
    write_minimal_files(run_dir)

    thresholds = mod.QCThresholds(
        ess_fail=150,
        ess_warn=180,
        overlap_fail=0.2,
        tilt_span=5.0,
        bulk_start_nm=2.0,
        bulk_min_points=3,
    )
    results = mod.analyze_simulation_runs(tmp_path, thresholds)
    out_dir = tmp_path / "qc_out"
    mod.write_reports(results, out_dir, mark_resim=True)

    csv_path = out_dir / "qc_report.csv"
    md_path = out_dir / "qc_report.md"
    yaml_path = out_dir / "qc_resimulate.yaml"
    assert csv_path.exists()
    assert md_path.exists()
    assert yaml_path.exists()
    flag_path = run_dir / "RESIMULATE.flag"
    assert flag_path.exists()
