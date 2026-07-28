#!/usr/bin/env python3
"""
Build a small genome subset + matching metadata file for fast local testing
of the annotate/enrichment/classify pixi tasks.

Truncates each 's__Faecalibacterium ...' group in data/groups.tsv down to
its first N genomes (sorted by genome id, for determinism) and keeps every
genome in all other (non-Faecalibacterium) groups untouched, since those
groups are already small.

Usage:
    python scripts/subset_genomes.py \
        --genome_dir data/genomes \
        --metadata data/groups.tsv \
        --output_genome_dir data/genomes_subset \
        --output_metadata data/groups_subset.tsv \
        --per_group 3
"""

import argparse
import os
import shutil


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--genome_dir', default='data/genomes')
    p.add_argument('--metadata', default='data/groups.tsv')
    p.add_argument('--output_genome_dir', default='data/genomes_subset')
    p.add_argument('--output_metadata', default='data/groups_subset.tsv')
    p.add_argument('--per_group', type=int, default=3,
                    help='Max genomes to keep per Faecalibacterium group (default: 3)')
    p.add_argument('--group_prefix', default='s__Faecalibacterium',
                    help='Only groups starting with this prefix are truncated (default: s__Faecalibacterium)')
    args = p.parse_args()

    groups = {}
    with open(args.metadata) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            genome, group = line.split('\t')
            groups.setdefault(group, []).append(genome)

    kept = []
    for group, genomes in groups.items():
        genomes = sorted(genomes)
        if group.startswith(args.group_prefix):
            genomes = genomes[:args.per_group]
        kept.extend((genome, group) for genome in genomes)

    os.makedirs(args.output_genome_dir, exist_ok=True)
    with open(args.output_metadata, 'w') as out:
        for genome, group in kept:
            out.write(f'{genome}\t{group}\n')

            src = os.path.join(args.genome_dir, f'{genome}.fna')
            dst = os.path.join(args.output_genome_dir, f'{genome}.fna')
            if not os.path.exists(src):
                raise FileNotFoundError(f'Genome file not found: {src}')
            if not os.path.exists(dst):
                shutil.copyfile(src, dst)

    print(f'Kept {len(kept)} genomes across {len(groups)} groups.')
    print(f'Genomes:  {args.output_genome_dir}/')
    print(f'Metadata: {args.output_metadata}')


if __name__ == '__main__':
    main()
