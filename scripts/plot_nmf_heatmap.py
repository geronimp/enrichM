#!/usr/bin/env python3
"""
Plot NMF genome scores as a heatmap annotated with group labels.

Usage:
    python scripts/plot_nmf_heatmap.py \
        --scores <nmf_scores.tsv> \
        --metadata <groups.tsv> \
        --output <output.png>

    pixi run nmf_heatmap --scores ... --metadata ... --output ...
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


def load_metadata(path):
    return pd.read_csv(path, sep='\t', header=None, names=['genome', 'group'])


def short_label(label):
    """Shorten GTDB-style labels for display: s__Faecalibacterium prausnitzii_F → F.prausnitzii_F"""
    label = label.replace('s__', '')
    parts = label.split()
    if len(parts) >= 2:
        return f'{parts[0][0]}.{parts[1]}'
    return label


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--scores', required=True, help='nmf_scores.tsv from enrichm enrichment')
    p.add_argument('--metadata', required=True, help='Genome-to-group TSV (no header)')
    p.add_argument('--output', default='nmf_heatmap.png', help='Output image path')
    p.add_argument('--figsize', nargs=2, type=float, default=[8, 12], metavar=('W', 'H'))
    p.add_argument('--dpi', type=int, default=150)
    args = p.parse_args()

    scores = pd.read_csv(args.scores, sep='\t', index_col=0)
    meta = load_metadata(args.metadata)
    meta = meta.set_index('genome')

    # Align metadata to scores, fill missing as 'Unknown'
    groups = meta.reindex(scores.index)['group'].fillna('Unknown')

    # Sort by group then by dominant component
    dominant = scores.idxmax(axis=1)
    sort_key = pd.DataFrame({'group': groups, 'dominant': dominant})
    order = sort_key.sort_values(['group', 'dominant']).index
    scores = scores.loc[order]
    groups = groups.loc[order]

    # Group colour palette
    unique_groups = groups.unique()
    palette = plt.cm.Set2(np.linspace(0, 0.8, len(unique_groups)))
    group_colours = {g: palette[i] for i, g in enumerate(unique_groups)}
    row_colours = groups.map(group_colours).values

    # --- Figure layout ---
    fig = plt.figure(figsize=args.figsize)
    # columns: group strip | heatmap | colourbar
    gs = fig.add_gridspec(1, 3, width_ratios=[0.04, 1, 0.04], wspace=0.02)

    ax_strip = fig.add_subplot(gs[0])
    ax_heat  = fig.add_subplot(gs[1])
    ax_cbar  = fig.add_subplot(gs[2])

    # Heatmap
    data = scores.values
    vmin, vmax = data.min(), data.max()
    im = ax_heat.imshow(data, aspect='auto', cmap='YlOrRd',
                        norm=Normalize(vmin=vmin, vmax=vmax),
                        interpolation='nearest')

    ax_heat.set_xticks(range(scores.shape[1]))
    ax_heat.set_xticklabels(scores.columns, rotation=45, ha='right', fontsize=9)
    ax_heat.set_yticks([])
    ax_heat.set_xlabel('NMF component', fontsize=10)
    ax_heat.set_title('NMF genome scores', fontsize=11, pad=8)

    # Group colour strip (left)
    strip = np.array([[matplotlib.colors.to_rgba(c)] for c in row_colours])
    ax_strip.imshow(strip, aspect='auto', interpolation='nearest')
    ax_strip.set_xticks([])
    ax_strip.set_yticks(range(len(order)))
    ax_strip.set_yticklabels(order, fontsize=5)
    ax_strip.yaxis.set_tick_params(length=0)
    ax_strip.set_ylabel('Genome', fontsize=9)

    # Colourbar for scores
    plt.colorbar(im, cax=ax_cbar, label='Score')

    # Legend for groups
    patches = [mpatches.Patch(color=group_colours[g], label=short_label(g))
               for g in unique_groups]
    ax_heat.legend(handles=patches, title='Group',
                   bbox_to_anchor=(1.18, 1), loc='upper left',
                   fontsize=8, title_fontsize=9, frameon=True)

    plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
    print(f'Saved: {args.output}')


if __name__ == '__main__':
    main()
