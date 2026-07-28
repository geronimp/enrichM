#!/usr/bin/env python3
"""
Run the enrichM annotate → enrichment pipeline in sequence.

Usage:
    pixi run pipeline --genome_dir <dir> --metadata <tsv> [options]

Options:
    --genome_dir    Directory of genome .fna files (required)
    --metadata      Genome-to-group TSV file (required)
    --output        Output directory prefix (default: pipeline_output)
    --threads       Number of threads for annotate (default: 8)
    --annotation    Annotation type: ko, ko_hmm, pfam, tigrfam, cazy, ec (default: ko)
    --decompose     Run NMF decomposition on enrichment results
    --tree          Path to Newick tree for phylogenetic correction
    --min_prevalence  Minimum genome prevalence fraction for annotations (default: 0.0)
    --force         Overwrite existing output directories
"""

import argparse
import subprocess
import sys
import os

def run(cmd):
    print(f'\n>>> {" ".join(cmd)}\n')
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(result.returncode)

def main():
    p = argparse.ArgumentParser(description='enrichM annotate + enrichment pipeline')
    p.add_argument('--genome_dir', required=True, help='Directory of genome .fna files')
    p.add_argument('--metadata', required=True, help='Genome-to-group TSV (no header)')
    p.add_argument('--output', default='pipeline_output', help='Output directory prefix (default: pipeline_output)')
    p.add_argument('--threads', type=int, default=8, help='Threads for annotate (default: 8)')
    p.add_argument('--annotation', default='ko',
                   choices=['ko', 'ko_hmm', 'pfam', 'tigrfam', 'cazy', 'ec'],
                   help='Annotation type to test (default: ko)')
    p.add_argument('--decompose', action='store_true', help='Run NMF decomposition')
    p.add_argument('--tree', default=None, help='Newick tree for phylogenetic correction')
    p.add_argument('--min_prevalence', type=float, default=0.0,
                   help='Minimum annotation prevalence fraction (default: 0.0)')
    p.add_argument('--force', action='store_true', help='Overwrite existing output directories')
    args = p.parse_args()

    annotate_dir = os.path.join(args.output, 'annotate')
    enrichment_dir = os.path.join(args.output, 'enrichment')

    enrichm = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'bin', 'enrichm')

    # --- Annotate ---
    annotate_cmd = [
        'python', enrichm, 'annotate',
        '--genome_directory', args.genome_dir,
        f'--{args.annotation}',
        '--threads', str(args.threads),
        '--output', annotate_dir,
    ]
    if args.force:
        annotate_cmd.append('--force')
    run(annotate_cmd)

    # --- Enrichment ---
    enrichment_cmd = [
        'python', enrichm, 'enrichment',
        '--annotate_output', annotate_dir,
        '--metadata', args.metadata,
        f'--{args.annotation}',
        '--min_prevalence', str(args.min_prevalence),
        '--output', enrichment_dir,
    ]
    if args.decompose:
        enrichment_cmd.append('--decompose')
    if args.tree:
        enrichment_cmd += ['--tree', args.tree]
    if args.force:
        enrichment_cmd.append('--force')
    run(enrichment_cmd)

    print(f'\nDone. Results in {enrichment_dir}/')

if __name__ == '__main__':
    main()
