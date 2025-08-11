#!/usr/bin/env python3
'''Extract features from RBC membrane simulation for SOLVIA_1'''

import numpy as np
import subprocess
import json

def extract_features():
    features = {}
    
    # Basic trajectory info
    cmd = ['gmx', 'check', '-f', 'md.xtc']
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Parse peptide-membrane distances
    if os.path.exists('peptide_membrane_dist.xvg'):
        distances = []
        with open('peptide_membrane_dist.xvg', 'r') as f:
            for line in f:
                if not line.startswith('#') and not line.startswith('@'):
                    parts = line.split()
                    if len(parts) >= 2:
                        distances.append(float(parts[1]))
        
        if distances:
            features['mean_distance'] = np.mean(distances)
            features['min_distance'] = np.min(distances)
            features['contact_frequency'] = sum(1 for d in distances if d < 0.5) / len(distances)
    
    # Save features
    with open('membrane_features.json', 'w') as f:
        json.dump(features, f, indent=2)
    
    print(f"Extracted features: {features}")
    return features

if __name__ == "__main__":
    extract_features()
