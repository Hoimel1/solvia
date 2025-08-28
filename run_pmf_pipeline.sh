#!/bin/bash
#
# SOLVIA PMF Master Pipeline
# Complete automated workflow for hemolytic toxicity prediction
#
# Usage: ./run_pmf_pipeline.sh <peptide_id> [options]
# Example: ./run_pmf_pipeline.sh 1

set -e  # Exit on error

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts/universal"
ANALYSIS_DIR="${SCRIPT_DIR}/scripts/analysis"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[PMF Pipeline]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Check arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 <peptide_id> [--skip-colabfold] [--replicates N] [--force]"
    echo "Example: $0 1"
    echo ""
    echo "Options:"
    echo "  --skip-colabfold  Skip structure prediction (use existing)"
    echo "  --replicates N    Number of replicates (default: 2)"
    echo "  --force           Force overwrite existing results"
    exit 1
fi

PEPTIDE_ID=$1
SKIP_COLABFOLD=0
N_REPLICATES=2
FORCE=0

# Parse additional arguments
shift
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-colabfold)
            SKIP_COLABFOLD=1
            shift
            ;;
        --replicates)
            N_REPLICATES=$2
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Setup paths
FASTA_FILE="${SCRIPT_DIR}/data/raw/fasta/SOLVIA_${PEPTIDE_ID}.fasta"
RUN_DIR="${SCRIPT_DIR}/simulations/solvia_${PEPTIDE_ID}_run_1"
MEMBRANE_DIR="${SCRIPT_DIR}/membrane_standard"

print_status "Starting PMF pipeline for SOLVIA_${PEPTIDE_ID}"
echo "=========================================="
echo "Configuration:"
echo "  FASTA: ${FASTA_FILE}"
echo "  Run directory: ${RUN_DIR}"
echo "  Replicates: ${N_REPLICATES}"
echo "  Skip ColabFold: ${SKIP_COLABFOLD}"
echo "=========================================="

# Check if FASTA exists
if [ ! -f "${FASTA_FILE}" ]; then
    print_error "FASTA file not found: ${FASTA_FILE}"
    exit 1
fi

# Check if run already exists
if [ -d "${RUN_DIR}" ] && [ ${FORCE} -eq 0 ]; then
    print_warning "Run directory already exists. Use --force to overwrite"
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# ============================================================
# STEP 1: Setup run directory
# ============================================================
print_status "Step 1: Setting up run directory"
python3 "${SCRIPTS_DIR}/01_setup_run.py" "${FASTA_FILE}"

# ============================================================
# STEP 2: Structure prediction (optional)
# ============================================================
if [ ${SKIP_COLABFOLD} -eq 0 ]; then
    print_status "Step 2: Running ColabFold structure prediction"
    # Check if ColabFold output already exists
    if [ -d "${RUN_DIR}/colabfold" ] && [ -f "${RUN_DIR}/colabfold/model_selection.yaml" ]; then
        print_warning "ColabFold results already exist, skipping"
    else
        python3 "${SCRIPTS_DIR}/02_run_colabfold.py" "${RUN_DIR}"
    fi
else
    print_status "Step 2: Skipping ColabFold (using existing structure)"
fi

# ============================================================
# STEP 3: Coarse-graining
# ============================================================
print_status "Step 3: Coarse-graining peptide structure"
if [ -f "${RUN_DIR}/coarse_grain/peptide_cg.gro" ]; then
    print_warning "Coarse-grained structure exists, skipping"
else
    python3 "${SCRIPTS_DIR}/03_coarse_grain.py" "${RUN_DIR}"
fi

# ============================================================
# STEP 4: Build standard membrane (if needed)
# ============================================================
print_status "Step 4: Preparing standard RBC membrane"
if [ ! -d "${MEMBRANE_DIR}" ]; then
    print_status "Building standard membrane template"
    python3 "${SCRIPTS_DIR}/04_build_membrane_standard.py" \
        --output-dir "${MEMBRANE_DIR}"
else
    print_status "Standard membrane already exists"
fi

# ============================================================
# STEP 5: Process each replicate
# ============================================================
for REP in $(seq 1 ${N_REPLICATES}); do
    echo ""
    print_status "Processing Replicate ${REP}/${N_REPLICATES}"
    echo "------------------------------------------"
    
    # 5a. Insert peptide with proper orientation
    print_status "Step 5a.${REP}: Inserting peptide (replicate ${REP})"
    if [ -d "${RUN_DIR}/pmf_system/replicate_${REP}" ]; then
        print_warning "PMF system for replicate ${REP} exists, skipping insertion"
    else
        python3 "${SCRIPTS_DIR}/05_insert_peptides_pmf.py" "${RUN_DIR}" \
            --replicate ${REP} \
            --peptide-id "SOLVIA_${PEPTIDE_ID}" \
            --posres-lipid-z
    fi
    
    # 5b. Equilibration
    print_status "Step 5b.${REP}: Running equilibration (replicate ${REP})"
    if [ -f "${RUN_DIR}/pmf_system/replicate_${REP}/equilibration/equilibration_summary.yaml" ]; then
        print_warning "Equilibration for replicate ${REP} complete, skipping"
    else
        python3 "${SCRIPTS_DIR}/06_equilibrate_pmf.py" "${RUN_DIR}" \
            --replicate ${REP}
    fi
    
    # 5c. PMF calculation
    print_status "Step 5c.${REP}: Running PMF umbrella sampling (replicate ${REP})"
    if [ -d "${RUN_DIR}/pmf/replicate_${REP}/windows" ]; then
        print_warning "PMF windows for replicate ${REP} exist, checking completeness..."
        # Could add check for window completeness here
    fi
    python3 "${SCRIPTS_DIR}/08_run_pmf_enhanced.py" "${RUN_DIR}" \
        --replicate ${REP}
    
    # 5d. Quality control
    print_status "Step 5d.${REP}: Running quality control (replicate ${REP})"
    python3 "${ANALYSIS_DIR}/pmf_qc_system.py" \
        "${RUN_DIR}/pmf/replicate_${REP}" \
        --suggest || true  # Don't exit on QC failure
    
    # 5e. MBAR analysis
    print_status "Step 5e.${REP}: Analyzing PMF with MBAR (replicate ${REP})"
    python3 "${ANALYSIS_DIR}/pmf_mbar_analysis.py" \
        "${RUN_DIR}/pmf/replicate_${REP}" \
        --method mbar \
        --bootstrap 1000
done

# ============================================================
# STEP 6: Aggregate results
# ============================================================
print_status "Step 6: Aggregating results from all replicates"

# Create results summary
RESULTS_DIR="${RUN_DIR}/final_results"
mkdir -p "${RESULTS_DIR}"

# Collect all replicate results
for REP in $(seq 1 ${N_REPLICATES}); do
    if [ -f "${RUN_DIR}/pmf/replicate_${REP}/pmf_analysis_results.yaml" ]; then
        cp "${RUN_DIR}/pmf/replicate_${REP}/pmf_analysis_results.yaml" \
           "${RESULTS_DIR}/replicate_${REP}_results.yaml"
    fi
done

# Create final summary
cat > "${RESULTS_DIR}/pipeline_summary.txt" << EOF
PMF Pipeline Summary for SOLVIA_${PEPTIDE_ID}
==============================================
Date: $(date)
Peptide ID: ${PEPTIDE_ID}
Replicates: ${N_REPLICATES}
Run Directory: ${RUN_DIR}

Results Files:
EOF

for REP in $(seq 1 ${N_REPLICATES}); do
    echo "  - Replicate ${REP}: replicate_${REP}_results.yaml" >> "${RESULTS_DIR}/pipeline_summary.txt"
done

# ============================================================
# STEP 7: Generate final report
# ============================================================
print_status "Step 7: Generating final report"

# Extract key features
echo "" >> "${RESULTS_DIR}/pipeline_summary.txt"
echo "Key PMF Features:" >> "${RESULTS_DIR}/pipeline_summary.txt"
echo "-----------------" >> "${RESULTS_DIR}/pipeline_summary.txt"

for REP in $(seq 1 ${N_REPLICATES}); do
    RESULTS_FILE="${RESULTS_DIR}/replicate_${REP}_results.yaml"
    if [ -f "${RESULTS_FILE}" ]; then
        echo "" >> "${RESULTS_DIR}/pipeline_summary.txt"
        echo "Replicate ${REP}:" >> "${RESULTS_DIR}/pipeline_summary.txt"
        
        # Extract features using Python
        python3 -c "
import yaml
with open('${RESULTS_FILE}', 'r') as f:
    data = yaml.safe_load(f)
    features = data.get('features', {})
    print(f\"  ΔG_ads: {features.get('delta_g_ads', 'N/A'):.2f} kJ/mol\" if features.get('delta_g_ads') is not None else '  ΔG_ads: N/A')
    print(f\"  ΔG_insert: {features.get('delta_g_insert', 'N/A'):.2f} kJ/mol\" if features.get('delta_g_insert') is not None else '  ΔG_insert: N/A')
    print(f\"  ΔG_barrier: {features.get('delta_g_barrier', 'N/A'):.2f} kJ/mol\" if features.get('delta_g_barrier') is not None else '  ΔG_barrier: N/A')
    
    # Check toxicity criteria
    if features.get('delta_g_ads') is not None:
        toxic = features.get('delta_g_ads') < -8.0
        if features.get('delta_g_insert') is not None:
            toxic = toxic and (features.get('delta_g_insert') <= -3.0 or 
                              (features.get('delta_g_barrier') is not None and features.get('delta_g_barrier') < 12.0))
        print(f\"  Toxicity prediction: {'TOXIC' if toxic else 'NON-TOXIC'}\")
" >> "${RESULTS_DIR}/pipeline_summary.txt"
    fi
done

# ============================================================
# FINAL STATUS
# ============================================================
echo ""
echo "=========================================="
print_status "PMF Pipeline Completed Successfully!"
echo "=========================================="
echo ""
echo "Results saved to: ${RESULTS_DIR}"
echo ""
echo "Key output files:"
echo "  - PMF profiles: ${RUN_DIR}/pmf/replicate_*/analysis_plots/pmf_profile.png"
echo "  - QC reports: ${RUN_DIR}/pmf/replicate_*/qc_full_report.yaml"
echo "  - Analysis results: ${RESULTS_DIR}/replicate_*_results.yaml"
echo "  - Summary: ${RESULTS_DIR}/pipeline_summary.txt"
echo ""

# Display summary
cat "${RESULTS_DIR}/pipeline_summary.txt"

exit 0
