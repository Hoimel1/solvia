"""
Visualization workflow for results
Based on evaluation and validation (Section 8)
"""

rule plot_roc_curve:
    input:
        predictions = DATA_DIR / "output" / "predictions" / "test_predictions.csv",
        metrics = DATA_DIR / "output" / "predictions" / "evaluation_metrics.json"
    output:
        roc_plot = DATA_DIR / "output" / "visualizations" / "roc_curve.png",
        pr_plot = DATA_DIR / "output" / "visualizations" / "pr_curve.png"
    log:
        "logs/visualization/roc_pr_curves.log"
    script:
        "../src/visualization/plot_curves.py"

rule plot_feature_importance:
    input:
        importance = DATA_DIR / "intermediate" / "ml" / "feature_importance.csv"
    output:
        plot = DATA_DIR / "output" / "visualizations" / "feature_importance.png"
    params:
        top_n = 20
    log:
        "logs/visualization/feature_importance.log"
    script:
        "../src/visualization/plot_feature_importance.py"

rule plot_prediction_distribution:
    input:
        predictions = DATA_DIR / "output" / "predictions" / "final_predictions.csv"
    output:
        dist_plot = DATA_DIR / "output" / "visualizations" / "prediction_distribution.png",
        calibration_plot = DATA_DIR / "output" / "visualizations" / "calibration_plot.png"
    log:
        "logs/visualization/prediction_distribution.log"
    script:
        "../src/visualization/plot_distributions.py"

rule generate_report:
    input:
        predictions = DATA_DIR / "output" / "predictions" / "final_predictions.csv",
        metrics = DATA_DIR / "output" / "predictions" / "evaluation_metrics.json",
        roc_plot = DATA_DIR / "output" / "visualizations" / "roc_curve.png",
        shap_summary = DATA_DIR / "output" / "visualizations" / "shap_summary.png",
        feature_importance = DATA_DIR / "output" / "visualizations" / "feature_importance.png"
    output:
        report = DATA_DIR / "output" / "report.html"
    params:
        config = config
    log:
        "logs/visualization/generate_report.log"
    script:
        "../src/visualization/generate_report.py"
