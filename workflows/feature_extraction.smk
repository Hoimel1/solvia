"""
Feature extraction workflow from MD trajectories
Based on feature sets (Section 4.2) and biophysical principles
"""

rule extract_features:
    input:
        traj = DATA_DIR / "intermediate" / "trajectories" / "{amp_id}" / "{occupancy}" / "rep{replicate}" / "prod.xtc",
        tpr = DATA_DIR / "intermediate" / "trajectories" / "{amp_id}" / "{occupancy}" / "rep{replicate}" / "prod.tpr",
        sequence = DATA_DIR / "input" / "sequences" / "{amp_id}.fasta"
    output:
        features = DATA_DIR / "intermediate" / "features" / "{amp_id}" / "{occupancy}" / "rep{replicate}" / "features.parquet"
    params:
        analysis_start = 400  # ns - analyze last 100 ns of 500 ns trajectory
    log:
        "logs/feature_extraction/{amp_id}_{occupancy}_rep{replicate}.log"
    script:
        "../src/feature_extraction/extract_md_features.py"

rule extract_sequence_features:
    input:
        sequence = DATA_DIR / "input" / "sequences" / "{amp_id}.fasta"
    output:
        features = DATA_DIR / "intermediate" / "features" / "{amp_id}" / "sequence_features.parquet"
    log:
        "logs/feature_extraction/{amp_id}_sequence.log"
    script:
        "../src/feature_extraction/extract_sequence_features.py"

rule aggregate_features:
    input:
        md_features = expand(
            DATA_DIR / "intermediate" / "features" / "{{amp_id}}" / "{occupancy}" / "rep{replicate}" / "features.parquet",
            occupancy=OCCUPANCIES,
            replicate=REPLICATES
        ),
        seq_features = DATA_DIR / "intermediate" / "features" / "{amp_id}" / "sequence_features.parquet"
    output:
        aggregated = DATA_DIR / "intermediate" / "features" / "{amp_id}" / "aggregated_features.parquet"
    log:
        "logs/feature_extraction/{amp_id}_aggregate.log"
    script:
        "../src/feature_extraction/aggregate_features.py"

rule combine_all_features:
    input:
        features = expand(
            DATA_DIR / "intermediate" / "features" / "{amp_id}" / "aggregated_features.parquet",
            amp_id=AMP_IDS
        ),
        labels = config["input"]["labels_file"]
    output:
        dataset = DATA_DIR / "intermediate" / "features" / "complete_dataset.parquet"
    log:
        "logs/feature_extraction/combine_all.log"
    run:
        import pandas as pd
        
        # Load all features
        all_features = []
        for feature_file in input.features:
            df = pd.read_parquet(feature_file)
            all_features.append(df)
        
        # Combine
        features_df = pd.concat(all_features, ignore_index=True)
        
        # Add labels
        labels_df = pd.read_csv(input.labels)
        final_df = features_df.merge(labels_df, on='amp_id')
        
        # Add binary classification
        final_df['is_toxic'] = (final_df['HC50'] <= config['ml']['toxicity_threshold']).astype(int)
        
        # Save
        final_df.to_parquet(output.dataset, index=False)

