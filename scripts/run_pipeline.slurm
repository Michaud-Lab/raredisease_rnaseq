#!/bin/bash
# =============================================================================
# run_pipeline.slurm - Full RNAseq rare-disease pipeline, run as one SLURM job
#
# Runs, in order: featureCounts -> OUTRIDER -> FRASER -> consensus -> ASE ->
# cp_and_cleanup.R
#
# Each step is a straight copy of the commands in its own <step>/<step>.slurm
# file (featureCounts/featureCounts.slurm, OUTRIDER/outrider.slurm, etc.),
# just run sequentially in one job instead of as separate submissions.
# Submit from anywhere with: sbatch run_pipeline.slurm
# =============================================================================

#SBATCH --job-name=raredisease_rnaseq_pipeline
#SBATCH --time=03:00:00
#SBATCH --mem=128G
#SBATCH --account=def-rallard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18

set -euo pipefail

PROJECT_ROOT="/project/def-rallard/COMMUN/raredisease_rnaseq"
SCRIPTS_DIR="$PROJECT_ROOT/scripts"
CONFIGS="$PROJECT_ROOT/data/input/configs.json"

cpu=16

# -----------------------------------------------------------------------------
# 0. Read every value used by any step out of configs.json in one jq call.
#    (featureCounts.fc_exons and consensus.fc_exons are different paths despite
#    the shared key name, so consensus' copy gets its own variable.)
# -----------------------------------------------------------------------------
read -r workdir rnasplice_bamdir gtf \
        fc_exons fc_genes gtf_out gtf_MANE ens_gene masterlog \
        fc_pergene fc_perexon \
        fraser_temp_bamdir fraser_bamdir fraser_chr_bamdir \
        fc_exons_tpm genome < <(
  jq -r '[
    .general.workdir, .general.rnasplice_bamdir, .general.gtf,
    .featureCounts.fc_exons, .featureCounts.fc_genes, .featureCounts.gtf_out,
      .featureCounts.gtf_MANE, .featureCounts.ens_gene, .featureCounts.masterlog,
    .OUTRIDER.fc_pergene, .OUTRIDER.fc_perexon,
    .FRASER.fraser_temp_bamdir, .FRASER.fraser_bamdir, .FRASER.fraser_chr_bamdir,
    .consensus.fc_exons_tpm, .consensus.genome
  ] | @tsv' "$CONFIGS"
)

# -----------------------------------------------------------------------------
# 1. featureCounts
# -----------------------------------------------------------------------------
echo "=== STEP 1/6: featureCounts ~~~ $(date) ==="
cd "$SCRIPTS_DIR/featureCounts"
module load subread/2.0.6
module load r

./featureCounts.sh $cpu $gtf $gtf_out $rnasplice_bamdir $workdir $gtf_MANE $ens_gene $masterlog $fc_exons $fc_genes
echo "=== DONE featureCounts ~~~ $(date) ==="

# -----------------------------------------------------------------------------
# 2. OUTRIDER
# -----------------------------------------------------------------------------
echo "=== STEP 2/6: OUTRIDER ~~~ $(date) ==="
cd "$SCRIPTS_DIR/OUTRIDER"

Rscript outrider.R $workdir $fc_pergene $fc_perexon $cpu
echo "=== DONE OUTRIDER ~~~ $(date) ==="

# -----------------------------------------------------------------------------
# 3. FRASER
# -----------------------------------------------------------------------------
echo "=== STEP 3/6: FRASER ~~~ $(date) ==="
cd "$SCRIPTS_DIR/FRASER"
module load samtools

Rscript fraser.R $workdir $rnasplice_bamdir $fraser_temp_bamdir $fraser_bamdir

# Genome-wide FRASER (fraser_gw.R) is left disabled here too, matching
# FRASER/fraser.slurm, where this loop is commented out:
#cpu=32
#for i in $(seq 1 22) X Y MT; do
#  Rscript fraser_gw.R $fraser_chr_bamdir $i $cpu $rnasplice_bamdir
#  echo 'Done chr:'$i
#done
echo "=== DONE FRASER ~~~ $(date) ==="

# -----------------------------------------------------------------------------
# 4. consensus
# -----------------------------------------------------------------------------
echo "=== STEP 4/6: consensus ~~~ $(date) ==="
cd "$SCRIPTS_DIR/consensus"
module load bedtools
module load bcftools

Rscript consensus.R $workdir $fc_exons_tpm $genome $gtf
echo "=== DONE consensus ~~~ $(date) ==="

# -----------------------------------------------------------------------------
# 5. ASE
# -----------------------------------------------------------------------------
echo "=== STEP 5/6: ASE ~~~ $(date) ==="
cd "$SCRIPTS_DIR/ASE"
module load gatk

./ASE.sh
echo "=== DONE ASE ~~~ $(date) ==="

# -----------------------------------------------------------------------------
# 6. cp_and_cleanup.R
# -----------------------------------------------------------------------------
echo "=== STEP 6/6: cp_and_cleanup ~~~ $(date) ==="
cd "$SCRIPTS_DIR"

Rscript cp_and_cleanup.R
echo "=== DONE cp_and_cleanup ~~~ $(date) ==="

echo "=== Pipeline finished ~~~ $(date) ==="
