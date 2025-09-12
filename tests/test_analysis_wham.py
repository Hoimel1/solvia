import os
from pathlib import Path
import numpy as np
import yaml


def make_synthetic_pmf_dir(tmp_path: Path, n_windows: int = 12) -> Path:
    # Create synthetic umbrella sampling dataset
    pmf_dir = tmp_path / "pmf_synth"
    pmf_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(1234)
    # Centers from -0.5 to 2.5 nm
    centers = np.linspace(-0.5, 2.5, n_windows)
    # Force constant and temperature must match analyzer defaults
    k = 900.0  # kJ/mol/nm^2
    T = 310.0

    windows = []
    for i, c in enumerate(centers):
        # Simple Gaussian around the center with sigma ~0.07 nm
        n = 1000
        t = np.arange(n, dtype=float)  # ps
        z = rng.normal(loc=c, scale=0.07, size=n)
        pullx_name = f"w{i:02d}_pullx.xvg"
        with open(pmf_dir / pullx_name, "w") as f:
            for ti, zi in zip(t, z):
                f.write(f"{ti:.3f} {zi:.5f}\n")
        windows.append({"center": float(c), "pullx": pullx_name})

    meta = {"windows": windows, "force_constant": float(k), "temperature": float(T)}
    with open(pmf_dir / "pmf_metadata.yaml", "w") as f:
        yaml.safe_dump(meta, f)
    return pmf_dir


def test_wham_default_runs_and_outputs(tmp_path):
    # Synthetic dataset
    pmf_dir = make_synthetic_pmf_dir(tmp_path)

    # Import analyzer module by path (no package required)
    import importlib.util
    mod_path = Path("scripts/analysis/pmf_mbar_analysis.py")
    spec = importlib.util.spec_from_file_location("pmf_mbar_analysis", mod_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)

    analyzer = module.PMFAnalyzer(pmf_dir)
    # Secure defaults: skip plots and bootstrap
    analyzer.config["no_plots"] = True
    analyzer.config["bootstrap"] = {"enabled": False, "n_bootstrap": 0}
    analyzer.config["method"] = "wham"

    results = analyzer.generate_analysis_report()
    assert results is not None

    # YAML written
    results_file = Path(pmf_dir) / "pmf_analysis_results.yaml"
    assert results_file.exists()
    with open(results_file, "r") as f:
        saved = yaml.safe_load(f)

    # Basic structure
    assert "pmf_data" in saved
    z = np.array(saved["pmf_data"]["z"], dtype=float)
    pmf = np.array(saved["pmf_data"]["pmf"], dtype=float)
    assert z.size == pmf.size and z.size > 8
    # Finite and not dominated by artificial walls
    assert np.isfinite(pmf).all()
    assert float(np.max(pmf)) < 1000.0


def test_cli_wham_default(tmp_path):
    pmf_dir = make_synthetic_pmf_dir(tmp_path)
    import subprocess, sys
    # Default CLI uses WHAM; disable plots/bootstraps
    cmd = [
        sys.executable,
        str(Path("scripts/analysis/pmf_mbar_analysis.py")),
        str(pmf_dir),
        "--no-plots",
        "--no-bootstrap",
    ]
    cp = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert cp.returncode == 0, f"CLI failed: {cp.stderr}\n{cp.stdout}"
    assert (Path(pmf_dir) / "pmf_analysis_results.yaml").exists()
