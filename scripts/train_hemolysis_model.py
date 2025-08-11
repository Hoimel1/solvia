#!/usr/bin/env python3
"""
Module 5: Train ML model for hemolysis prediction.

Uses XGBoost with SHAP for interpretability.
Includes cross-validation, hyperparameter tuning, and feature importance analysis.
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path
import logging
import json
import pickle
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.metrics import confusion_matrix, classification_report
import xgboost as xgb
import shap
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class HemolysisPredictor:
    def __init__(self, feature_file, label_file=None):
        self.feature_file = Path(feature_file)
        self.label_file = Path(label_file) if label_file else None
        self.model = None
        self.scaler = StandardScaler()
        self.feature_names = None
        
    def load_data(self):
        """Load features and labels."""
        logger.info("Loading data...")
        
        # Load features
        if self.feature_file.suffix == '.csv':
            self.df_features = pd.read_csv(self.feature_file)
        else:
            # Load from individual JSON files
            feature_dir = self.feature_file.parent
            all_features = []
            
            for json_file in feature_dir.glob("*_features.json"):
                with open(json_file, 'r') as f:
                    data = json.load(f)
                    feature_dict = {'peptide_id': data['peptide_id']}
                    feature_dict.update(data['features'])
                    all_features.append(feature_dict)
                    
            self.df_features = pd.DataFrame(all_features)
        
        # Load labels
        if self.label_file and self.label_file.exists():
            self.df_labels = pd.read_csv(self.label_file)
            
            # Merge features with labels
            self.df = pd.merge(
                self.df_features, 
                self.df_labels[['ID', 'Hemolysis']], 
                left_on='peptide_id', 
                right_on='ID',
                how='inner'
            )
            
            # Binary classification: Hemolytic (1) vs Non-hemolytic (0)
            self.df['hemolytic'] = (self.df['Hemolysis'] > 10).astype(int)
            
        else:
            logger.warning("No label file provided. Generating synthetic labels for demo.")
            # Generate synthetic labels based on features
            self.df = self.df_features.copy()
            
            # Simple rule-based synthetic labels
            self.df['hemolytic'] = (
                (self.df.get('membrane_contacts_mean', 50) > 60) & 
                (self.df.get('membrane_perturbation', 0.5) > 0.6)
            ).astype(int)
            
            # Add some noise
            np.random.seed(42)
            noise_idx = np.random.choice(len(self.df), size=int(0.1 * len(self.df)), replace=False)
            self.df.loc[noise_idx, 'hemolytic'] = 1 - self.df.loc[noise_idx, 'hemolytic']
            
        logger.info(f"Loaded {len(self.df)} samples")
        logger.info(f"Class distribution: {self.df['hemolytic'].value_counts().to_dict()}")
        
        # Prepare features
        self.feature_names = [col for col in self.df.columns 
                             if col not in ['peptide_id', 'ID', 'Hemolysis', 'hemolytic']]
        
        self.X = self.df[self.feature_names].values
        self.y = self.df['hemolytic'].values
        
        # Handle missing values
        self.X = np.nan_to_num(self.X, nan=0.0)
        
    def split_data(self, test_size=0.2, random_state=42):
        """Split data into train and test sets."""
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
            self.X, self.y, test_size=test_size, random_state=random_state, stratify=self.y
        )
        
        # Scale features
        self.X_train_scaled = self.scaler.fit_transform(self.X_train)
        self.X_test_scaled = self.scaler.transform(self.X_test)
        
        logger.info(f"Train set: {len(self.X_train)} samples")
        logger.info(f"Test set: {len(self.X_test)} samples")
        
    def train_model(self, use_grid_search=True):
        """Train XGBoost model with optional hyperparameter tuning."""
        logger.info("Training XGBoost model...")
        
        if use_grid_search and len(self.X_train) > 50:
            # Hyperparameter tuning
            param_grid = {
                'n_estimators': [50, 100, 200],
                'max_depth': [3, 5, 7],
                'learning_rate': [0.01, 0.1, 0.3],
                'subsample': [0.8, 1.0],
                'colsample_bytree': [0.8, 1.0]
            }
            
            xgb_model = xgb.XGBClassifier(
                objective='binary:logistic',
                random_state=42,
                use_label_encoder=False,
                eval_metric='logloss'
            )
            
            grid_search = GridSearchCV(
                xgb_model, 
                param_grid, 
                cv=5, 
                scoring='roc_auc',
                n_jobs=-1,
                verbose=1
            )
            
            grid_search.fit(self.X_train_scaled, self.y_train)
            self.model = grid_search.best_estimator_
            
            logger.info(f"Best parameters: {grid_search.best_params_}")
            logger.info(f"Best CV score: {grid_search.best_score_:.3f}")
            
        else:
            # Simple training
            self.model = xgb.XGBClassifier(
                n_estimators=100,
                max_depth=5,
                learning_rate=0.1,
                objective='binary:logistic',
                random_state=42,
                use_label_encoder=False,
                eval_metric='logloss'
            )
            
            self.model.fit(
                self.X_train_scaled, 
                self.y_train,
                eval_set=[(self.X_test_scaled, self.y_test)],
                verbose=False
            )
            
    def evaluate_model(self):
        """Evaluate model performance."""
        logger.info("Evaluating model...")
        
        # Predictions
        self.y_pred = self.model.predict(self.X_test_scaled)
        self.y_pred_proba = self.model.predict_proba(self.X_test_scaled)[:, 1]
        
        # Metrics
        metrics = {
            'accuracy': accuracy_score(self.y_test, self.y_pred),
            'precision': precision_score(self.y_test, self.y_pred),
            'recall': recall_score(self.y_test, self.y_pred),
            'f1': f1_score(self.y_test, self.y_pred),
            'auc': roc_auc_score(self.y_test, self.y_pred_proba)
        }
        
        logger.info("Model performance:")
        for metric, value in metrics.items():
            logger.info(f"  {metric}: {value:.3f}")
            
        # Confusion matrix
        cm = confusion_matrix(self.y_test, self.y_pred)
        logger.info(f"\nConfusion Matrix:\n{cm}")
        
        # Classification report
        logger.info("\nClassification Report:")
        logger.info(classification_report(self.y_test, self.y_pred, 
                                        target_names=['Non-hemolytic', 'Hemolytic']))
        
        # Cross-validation
        cv_scores = cross_val_score(self.model, self.X_train_scaled, self.y_train, 
                                   cv=5, scoring='roc_auc')
        logger.info(f"\n5-fold CV AUC: {cv_scores.mean():.3f} (+/- {cv_scores.std() * 2:.3f})")
        
        return metrics
        
    def analyze_features(self, output_dir='results'):
        """Analyze feature importance using SHAP."""
        logger.info("Analyzing feature importance...")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # XGBoost feature importance
        importance = self.model.feature_importances_
        indices = np.argsort(importance)[::-1]
        
        # Plot feature importance
        plt.figure(figsize=(10, 8))
        top_n = min(20, len(self.feature_names))
        plt.barh(range(top_n), importance[indices[:top_n]])
        plt.yticks(range(top_n), [self.feature_names[i] for i in indices[:top_n]])
        plt.xlabel('Feature Importance')
        plt.title('Top Features by XGBoost Importance')
        plt.tight_layout()
        plt.savefig(output_dir / 'feature_importance.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # SHAP analysis
        try:
            explainer = shap.TreeExplainer(self.model)
            shap_values = explainer.shap_values(self.X_test_scaled)
            
            # Summary plot
            plt.figure(figsize=(10, 8))
            shap.summary_plot(shap_values, self.X_test_scaled, 
                            feature_names=self.feature_names, 
                            show=False)
            plt.tight_layout()
            plt.savefig(output_dir / 'shap_summary.png', dpi=300, bbox_inches='tight')
            plt.close()
            
            # Feature importance based on SHAP
            shap_importance = np.abs(shap_values).mean(axis=0)
            
            # Save feature importance
            importance_df = pd.DataFrame({
                'feature': self.feature_names,
                'xgb_importance': importance,
                'shap_importance': shap_importance
            }).sort_values('shap_importance', ascending=False)
            
            importance_df.to_csv(output_dir / 'feature_importance.csv', index=False)
            
            logger.info("\nTop 10 features by SHAP importance:")
            for _, row in importance_df.head(10).iterrows():
                logger.info(f"  {row['feature']}: {row['shap_importance']:.3f}")
                
        except Exception as e:
            logger.warning(f"SHAP analysis failed: {e}")
            
    def save_model(self, output_dir='models'):
        """Save trained model and scaler."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save model
        model_path = output_dir / 'hemolysis_model.pkl'
        with open(model_path, 'wb') as f:
            pickle.dump(self.model, f)
            
        # Save scaler
        scaler_path = output_dir / 'feature_scaler.pkl'
        with open(scaler_path, 'wb') as f:
            pickle.dump(self.scaler, f)
            
        # Save feature names
        feature_path = output_dir / 'feature_names.json'
        with open(feature_path, 'w') as f:
            json.dump(self.feature_names, f)
            
        # Save model info
        info = {
            'model_type': 'XGBoost',
            'n_features': len(self.feature_names),
            'feature_names': self.feature_names,
            'performance': self.evaluate_model()
        }
        
        info_path = output_dir / 'model_info.json'
        with open(info_path, 'w') as f:
            json.dump(info, f, indent=2)
            
        logger.info(f"Model saved to {output_dir}")
        
    def predict_new_peptides(self, feature_file):
        """Predict hemolysis for new peptides."""
        # Load features
        if feature_file.endswith('.json'):
            with open(feature_file, 'r') as f:
                data = json.load(f)
                features = data['features']
        else:
            df = pd.read_csv(feature_file)
            features = df[self.feature_names].iloc[0].to_dict()
            
        # Prepare features
        X = np.array([[features.get(fn, 0) for fn in self.feature_names]])
        X_scaled = self.scaler.transform(X)
        
        # Predict
        prob = self.model.predict_proba(X_scaled)[0, 1]
        pred = self.model.predict(X_scaled)[0]
        
        return {
            'prediction': 'Hemolytic' if pred == 1 else 'Non-hemolytic',
            'probability': prob,
            'confidence': max(prob, 1 - prob)
        }


def create_synthetic_labels(feature_dir, output_file):
    """Create synthetic labels for demo purposes."""
    logger.info("Creating synthetic labels...")
    
    # Load peptide IDs from features
    peptide_ids = []
    for json_file in Path(feature_dir).glob("*_features.json"):
        peptide_id = json_file.stem.replace('_features', '')
        peptide_ids.append(peptide_id)
    
    # Generate synthetic hemolysis values
    np.random.seed(42)
    labels = []
    
    for peptide_id in peptide_ids:
        # Extract peptide number
        try:
            num = int(peptide_id.split('_')[1])
        except:
            num = np.random.randint(1, 2000)
            
        # Create pattern: lower numbers tend to be less hemolytic
        base_hemolysis = (num / 2000) * 50 + np.random.normal(0, 10)
        hemolysis = max(0, min(100, base_hemolysis))
        
        labels.append({
            'ID': peptide_id,
            'Hemolysis': hemolysis
        })
    
    df_labels = pd.DataFrame(labels)
    df_labels.to_csv(output_file, index=False)
    
    logger.info(f"Created synthetic labels for {len(labels)} peptides")
    logger.info(f"Hemolysis range: {df_labels['Hemolysis'].min():.1f} - {df_labels['Hemolysis'].max():.1f}")
    
    return df_labels


def main():
    """Run the training pipeline."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Train hemolysis prediction model")
    parser.add_argument('--features', default='data/features/all_features.csv', 
                       help='Feature file or directory')
    parser.add_argument('--labels', default='data/processed/hemolysis_labels.csv',
                       help='Label file')
    parser.add_argument('--create-labels', action='store_true',
                       help='Create synthetic labels for demo')
    parser.add_argument('--output-dir', default='results',
                       help='Output directory')
    parser.add_argument('--model-dir', default='models',
                       help='Model save directory')
    
    args = parser.parse_args()
    
    # Create synthetic labels if needed
    if args.create_labels or not Path(args.labels).exists():
        if Path(args.features).is_dir():
            feature_dir = args.features
        else:
            feature_dir = Path(args.features).parent
            
        args.labels = 'data/processed/synthetic_labels.csv'
        create_synthetic_labels(feature_dir, args.labels)
    
    # Initialize predictor
    predictor = HemolysisPredictor(args.features, args.labels)
    
    # Load and prepare data
    predictor.load_data()
    predictor.split_data()
    
    # Train model
    predictor.train_model(use_grid_search=False)  # Fast training for demo
    
    # Evaluate
    metrics = predictor.evaluate_model()
    
    # Analyze features
    predictor.analyze_features(args.output_dir)
    
    # Save model
    predictor.save_model(args.model_dir)
    
    logger.info("\nTraining complete!")
    logger.info(f"Results saved to {args.output_dir}")
    logger.info(f"Model saved to {args.model_dir}")


if __name__ == "__main__":
    main()
