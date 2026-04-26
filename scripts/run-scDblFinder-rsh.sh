#!/bin/bash
#SBATCH --job-name=scDblFinder_TIL
#SBATCH --nodes=1
#SBATCH --time=0-24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=1
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-user=jlwalker@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --account=rbao
#SBATCH --mem=800G
#SBATCH --output=/ix/rbao/Projects/HCC-CBS-253-Hillman-LKarapetyan-scRNAseq-revise/results/jwork/logs/scDblFinder_%j.out
#SBATCH --error=/ix/rbao/Projects/HCC-CBS-253-Hillman-LKarapetyan-scRNAseq-revise/results/jwork/logs/scDblFinder_%j.err

# ── Environment ──────────────────────────────────────────────────────────
module load gcc/12.2.0
module load r/4.5.0  # adjust to your R version

# ── Create log dir if needed ─────────────────────────────────────────────
mkdir -p /ix/rbao/Projects/HCC-CBS-253-Hillman-LKarapetyan-scRNAseq-revise/results/jwork/logs

# ── Run ──────────────────────────────────────────────────────────────────
cd /ix/rbao/Projects/HCC-CBS-253-Hillman-LKarapetyan-scRNAseq-revise/results/jwork

Rscript run_scDblFinder.R