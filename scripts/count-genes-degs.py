count-genes-degs.py

import pandas as pd
import re

# ── CONFIG ────────────────────────────────────────────────────────────────────
input_csv  = 'all-pfs-PBMC-combined-data.csv'   # swap for TIL file as needed
out_annotated   = 'all_pfs_PBMC_DEG_annotated.csv'
out_gene_counts = 'gene_counts_summary_PBMC.csv'
# ─────────────────────────────────────────────────────────────────────────────

df = pd.read_csv(input_csv)

# --- 1. Extract cell type from source_file path ---
def extract_cell_type(src):
    m = re.search(r'NoPrg6M\.vs\.Prg6M\.(.+?)\.DEG', src)
    if m:
        return m.group(1).replace('.', ' ')
    return None

df['cell.type']      = df['source_file'].apply(extract_cell_type)
df['cell.contrasts'] = df['contrast'] + '.' + df['cell.type'].str.replace(' ', '.')

# --- 2. Add summary column based on fold-change direction ---
# Contrast is NoPrg6M vs Prg6M, so:
#   avg_log2FC > 0  →  up in non-progressor
#   avg_log2FC < 0  →  up in progressor
df['summary'] = df['avg_log2FC'].apply(
    lambda x: 'PFS.6month.up.in.non.progressor' if x > 0
              else 'PFS.6month.up.in.progressor'
)

# --- 3. Select final columns for annotated DEG table ---
all_de = df[[
    'Gene', 'p_val_adj', 'avg_log2FC', 'summary',
    'pct.1', 'pct.2', 'contrast', 'cell.type',
    'cell.contrasts', 'diffexpressed'
]].copy()

# --- 4. Gene count table: how many times each gene appears ---
gene_counts = (
    all_de
    .groupby('Gene')
    .size()
    .reset_index(name='counts')
    .sort_values('counts', ascending=False)
)

# --- 5. Directionality breakdown per gene ---
dir_counts = (
    all_de
    .groupby(['Gene', 'summary'])
    .size()
    .reset_index(name='n')
    .pivot(index='Gene', columns='summary', values='n')
    .fillna(0)
    .astype(int)
    .reset_index()
)

# --- 6. Merge counts + directionality ---
gene_summary = gene_counts.merge(dir_counts, on='Gene', how='left')

# --- 7. Save outputs ---
all_de.to_csv(out_annotated, index=False)
gene_summary.to_csv(out_gene_counts, index=False)

print(f"Annotated DEG table:  {all_de.shape}  →  {out_annotated}")
print(f"Gene counts summary:  {gene_summary.shape}  →  {out_gene_counts}")
print(gene_summary.head(20).to_string(index=False))

#for TIL
input_csv       = 'all-pfs-TIL-combined-data.csv'
out_annotated   = 'all_pfs_TIL_DEG_annotated.csv'
out_gene_counts = 'gene_counts_summary_TIL.csv'

#regex pattern 'NoPrg6M\.vs\.Prg6M\.(.+?)\.DEG'

#recurring deg plots:
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib import cm
import re

# --- Filter function: removes uncharacterized / artifact genes ---
def is_artifact(gene):
    histone_prefixes = ('HIST', 'H1F', 'H2AF','H3F', 'H4F', 'H1-', 'H2A', 'H2B', 'H3-', 'H4-', 'H1.')
    return (
        gene.startswith('ENSG') or                      # unannotated Ensembl IDs
        gene.startswith('MT-') or                       # mitochondrial genes
        gene.startswith('MTRN') or                      # mitochondrial RNA genes
        any(gene.startswith(p) for p in histone_prefixes) or  # histone genes
        bool(re.match(r'^AC\d+', gene)) or              # uncharacterized AC loci
        gene.startswith('LINC') or                      # long intergenic non-coding RNAs
        bool(re.search(r'-AS\d+$', gene))               # antisense transcripts
    )

# --- Plot function ---
def make_plot(label, annotated_path, contrast_levels, out_png, out_pdf):
    all_de = pd.read_csv(annotated_path)

    # Apply artifact filter
    all_de = all_de[~all_de['Gene'].apply(is_artifact)].copy()

    # Top 50 genes by count, ascending order for coord_flip effect
    gene_order = (
        all_de.groupby('Gene')
        .size()
        .reset_index(name='counts')
        .sort_values('counts', ascending=False)
        .head(50)
        .sort_values('counts', ascending=True)
    )
    top_genes = gene_order['Gene'].tolist()

    # Filter to top genes and set categorical order
    temp = all_de[all_de['Gene'].isin(top_genes)].copy()
    temp['Gene'] = pd.Categorical(temp['Gene'], categories=top_genes, ordered=True)

    # Colors: Blues for Baseline contrasts, YlOrRd for Week4 contrasts
    blues = [cm.Blues(x) for x in np.linspace(0.3, 0.95, 9)]
    ylrd  = [cm.YlOrRd(x) for x in np.linspace(0.2, 0.95, 9)]
    cols  = blues + ylrd

    # Pivot to stacked bar format
    pivot = (
        temp.groupby(['Gene', 'cell.contrasts'])
        .size()
        .reset_index(name='n')
        .pivot(index='Gene', columns='cell.contrasts', values='n')
        .fillna(0)
        .reindex(index=top_genes)
        .reindex(columns=contrast_levels, fill_value=0)
    )

    # Plot
    fig, ax = plt.subplots(figsize=(9, 14))
    bottom = np.zeros(len(top_genes))
    y_pos  = np.arange(len(top_genes))

    for contrast, color in zip(contrast_levels, cols):
        if contrast in pivot.columns:
            vals = pivot[contrast].values
            ax.barh(y_pos, vals, left=bottom, color=color, height=0.7)
            bottom += vals

    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_genes, fontsize=8)
    ax.set_xlabel('Number of times found', fontsize=11)
    ax.set_xlim(0, max(bottom) + 1)
    ax.set_title(f'Top 50 Recurring DEGs — {label}\n(ENS, MT, Histone, AC, LINC, -AS genes excluded)',
                 fontsize=12, fontweight='bold')

    # Legend with shortened contrast labels
    short_labels = [
        c.replace(f'{label}.PFS.6month.NoPrg6M.vs.Prg6M.', '').replace('.', ' ')
        for c in contrast_levels
    ]
    patches = [mpatches.Patch(color=cols[i], label=short_labels[i])
               for i in range(len(contrast_levels))]
    ax.legend(handles=patches, bbox_to_anchor=(1.01, 1), loc='upper left',
              fontsize=7, frameon=False, title='Contrast', title_fontsize=8)

    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.savefig(out_pdf, bbox_inches='tight')
    plt.close()
    print(f'{label} saved.')


# --- Contrast levels ---
pbmc_contrasts = [
    'Baseline.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell',
    'Baseline.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages',
    'Baseline.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg',
    'Baseline.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell',
    'Baseline.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages',
    'Baseline.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg',
    'Baseline.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell',
    'Baseline.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages',
    'Baseline.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg',
    'Week4.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell',
    'Week4.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages',
    'Week4.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg',
    'Week4.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell',
    'Week4.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages',
    'Week4.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg',
    'Week4.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell',
    'Week4.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages',
    'Week4.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg',
]

# TIL contrasts are identical but with 'TIL' instead of 'PBMC'
til_contrasts = [c.replace('PBMC', 'TIL') for c in pbmc_contrasts]

# --- Run for both datasets ---
make_plot('PBMC',
          'all_pfs_PBMC_DEG_annotated.csv',
          pbmc_contrasts,
          'recurring_DEG_PBMC_final2.png',
          'recurring_DEG_PBMC_final2.pdf')

make_plot('TIL',
          'all_pfs_TIL_DEG_annotated.csv',
          til_contrasts,
          'recurring_DEG_TIL_final2.png',
          'recurring_DEG_TIL_final2.pdf')