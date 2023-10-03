"""
To benchmark on a set of proteins:
First, you will need to run the combind pipeline up to the "combind featurize"
step for all of the PDB ligands, i.e. those for which there are experimental
structures.

Then, use the "compute" function to compute the distributions for each
individual protein in the dataset. Here the feature_paths argument should
lead to the output from "combind featurize" runs.

Finally, create merged distributions for each protein that you want to benchmark
on using the "merge" function. Here the to_merge argument should include the
paths to all of the distributions except the protein that you are benchmarking
on.
"""

import sys
sys.path.append('/home/users/psuriana/combind-orig')

from score import statistics
from score.density_estimate import DensityEstimate
from features.features import Features
import os
import click
import pandas as pd
import numpy as np
from tqdm import tqdm
from glob import glob
from pathlib import Path


@click.group()
def main():
    pass


# JOE's COMBIND DATASET (from /oak/stanford/groups/rondror/users/jpaggi/combind_paper_dataset)
# python -m src.combind.refit_stats compute /scratch/groups/rondror/psuriana/combind_paper_dataset/stats/sp_es4 /scratch/groups/rondror/psuriana/combind_paper_dataset/*/stats_features/sp_es4 --protein-name-position -3

# python -m src.combind.refit_stats compute /scratch/users/psuriana/combind_paper_systems/stats/e3nn-run-20220422_101444-zybjdx87-best_val-max_helper_20 /scratch/users/psuriana/combind_paper_systems/data/*/stats_features/e3nn-run-20220422_101444-zybjdx87-best_val-max_helper_20 --protein-name-position -3

# TRAIN WITH TRAIN AND VAL SET
# python -m src.combind.refit_stats compute /scratch/users/psuriana/combind_paper_systems/stats/e3nn-run-20220509_115644-29cuk25e-max_helper_20 /scratch/users/psuriana/combind_paper_systems/data/*/stats_features/e3nn-run-20220509_115644-29cuk25e-max_helper_20 --protein-name-position -3

# TRAIN WITH ALL TRAIN SET
# python -m src.combind.refit_stats compute /scratch/users/psuriana/combind_paper_systems/stats/e3nn-run-20220508_160135-2kxspygq-max_helper_20 /scratch/users/psuriana/combind_paper_systems/data/*/stats_features/e3nn-run-20220508_160135-2kxspygq-max_helper_20 --protein-name-position -3

# GLIDE ORIGINAL
# python -m src.combind.refit_stats compute /scratch/users/psuriana/combind_paper_systems/stats/sp_es4-max_poses_300 /scratch/users/psuriana/combind_paper_systems/data/*/stats_features/sp_es4-max_poses_300-max_helper_20 --protein-name-position -3
#
# python -m src.combind.refit_stats compute /scratch/users/psuriana/combind_paper_systems/stats/e3nn_less_strict-run-20220422_171334-99eyoyz8-max_helper_20 /scratch/users/psuriana/combind_paper_systems/data/*/stats_features/e3nn_less_strict-run-20220422_171334-99eyoyz8-max_helper_20 --protein-name-position -3
@main.command()
@click.argument('stats_root')
@click.argument('feature_paths', nargs=-1)
@click.option('--protein-name-position', default=-2)
@click.option('--feature-names', default='hbond,saltbridge,contact,shape,mcss')
def compute(stats_root, protein_name_position, feature_paths, feature_names):
    feature_names = feature_names.split(',')
    all_proteins = set()
    for path in tqdm(feature_paths):
        print(path)
        protein = path.split('/')[protein_name_position]
        assert protein not in all_proteins, f'Found duplicate {protein}'
        all_proteins.add(protein)

        output_dir = os.path.join(stats_root, protein)
        os.makedirs(output_dir, exist_ok=True)
        print(f'Writing output to {output_dir}')

        features = Features(path)
        features.load_features()

        # Pairs that are between the same ligand
        names = features.raw['name1']
        self_self = names.reshape(-1, 1) == names.reshape(1, -1)

        # Avoid double counting pose pairs
        idx = np.array(range(len(features.raw['name1'])))
        upper_triangular = idx.reshape(-1, 1) < idx.reshape(1, -1)

        # Near-native pose pairs
        native = features.raw['rmsd1'] <= 2.0
        native = native.reshape(-1, 1) * native.reshape(1, -1)

        for feature in feature_names:
            if not feature in features.raw:
                continue
            print(f'Density estimate for {feature}')

            # Mask off value that is infinite (MCSS failed to compute)
            not_infinite = (features.raw[feature] != float('inf'))

            nat_vals = features.raw[feature][~self_self & upper_triangular & native & not_infinite]
            ref_vals = features.raw[feature][~self_self & upper_triangular & not_infinite]
            assert np.all(nat_vals != float('inf'))
            assert np.all(ref_vals != float('inf'))
            if feature == 'mcss':
                sd = 0.03*6
                #domain = (0, 4)
                # Combind used 6 as cutoff but I was unlucky and have a slight bump around 5
                # for native reference which skew the statistics, that's why we changed the
                # cutoff to 4 instead of 6
                domain = (0, 6)
            else:
                sd = 0.03
                domain = (0, 1)
            native_density = os.path.join(output_dir, f'native_{feature}.de')
            if (not os.path.exists(native_density)) or (feature == 'mcss'):
                nat = DensityEstimate(domain=domain, sd=sd).fit(nat_vals)
                print(f'Writing native density to {native_density}')
                nat.write(native_density)

            reference_density = os.path.join(output_dir, f'reference_{feature}.de')
            if (not os.path.exists(reference_density)) or (feature == 'mcss'):
                ref = DensityEstimate(domain=domain, sd=sd).fit(ref_vals)
                print(f'Writing reference density to {reference_density}')
                ref.write(reference_density)


# JOE's COMBIND DATASET (from /oak/stanford/groups/rondror/users/jpaggi/combind_paper_dataset)
# python -m src.combind.refit_stats merge /scratch/groups/rondror/psuriana/combind_paper_dataset/merged_stats/sp_es4/all /scratch/groups/rondror/psuriana/combind_paper_dataset/stats/sp_es4/*


# python -m src.combind.refit_stats merge /scratch/users/psuriana/combind_paper_systems/merged_stats/e3nn-run-20220422_101444-zybjdx87-best_val-max_helper_20/all /scratch/users/psuriana/combind_paper_systems/stats/e3nn-run-20220422_101444-zybjdx87-best_val-max_helper_20/*

# python -m src.combind.refit_stats merge /scratch/users/psuriana/combind_paper_systems/merged_stats/e3nn_less_strict-run-20220422_171334-99eyoyz8-max_helper_20/all /scratch/users/psuriana/combind_paper_systems/stats/e3nn_less_strict-run-20220422_171334-99eyoyz8-max_helper_20/*

# TRAIN WITH TRAIN AND VAL SET
# python -m src.combind.refit_stats merge /scratch/users/psuriana/combind_paper_systems/merged_stats/e3nn-run-20220509_115644-29cuk25e-max_helper_20/all /scratch/users/psuriana/combind_paper_systems/stats/e3nn-run-20220509_115644-29cuk25e-max_helper_20/*


# TRAIN WITH ALL TRAIN SET
# python -m src.combind.refit_stats merge /scratch/users/psuriana/combind_paper_systems/merged_stats/e3nn-run-20220508_160135-2kxspygq-max_helper_20/all /scratch/users/psuriana/combind_paper_systems/stats/e3nn-run-20220508_160135-2kxspygq-max_helper_20/*
## EXCLUDE TRANSPORTER
# python -m src.combind.refit_stats merge /scratch/users/psuriana/combind_paper_systems/merged_stats/e3nn-run-20220508_160135-2kxspygq-max_helper_20-exclude_transporters/all /scratch/users/psuriana/combind_paper_systems/stats/e3nn-run-20220508_160135-2kxspygq-max_helper_20/* -exclude SLC6A4 -exclude GLUT1 -exclude DAT

# GLIDE ORIGINAL
# python -m src.combind.refit_stats merge /scratch/users/psuriana/combind_paper_systems/merged_stats/sp_es4-max_poses_300/all /scratch/users/psuriana/combind_paper_systems/stats/sp_es4-max_poses_300/*
## EXCLUDE TRANSPORTER
# python -m src.combind.refit_stats merge /scratch/users/psuriana/combind_paper_systems/merged_stats/sp_es4-max_poses_300-exclude_transporters/all /scratch/users/psuriana/combind_paper_systems/stats/sp_es4-max_poses_300/* -exclude SLC6A4 -exclude GLUT1 -exclude DAT
@main.command()
@click.argument('merged')
@click.argument('to_merge', nargs=-1)
@click.option('--feature-names', default='hbond,saltbridge,contact,shape,mcss')
@click.option('--proteins_to_exclude', '-exclude', type=str, multiple=True,
              default=[] # Transporter: ['SLC6A4', 'GLUT1', 'DAT']
              )
def merge(merged, to_merge, feature_names, proteins_to_exclude):
    return merge_stats(merged, to_merge, feature_names, proteins_to_exclude)


def merge_stats(merged, to_merge, feature_names, proteins_to_exclude=[]):
    print(f'Find {len(to_merge)} files to process')
    if len(proteins_to_exclude) > 0:
        print(f'Exclude protein {proteins_to_exclude}')
        to_merge = [f for f in to_merge if Path(f).name not in proteins_to_exclude]
        print(f'After filtering find {len(to_merge)} files to process')
    os.makedirs(merged, exist_ok=True)
    feature_names = feature_names.split(',')
    for feature in feature_names:
        print(f'\nFeature {feature}')
        nat_des, ref_des = [], []
        for stats_root in to_merge:
            print(f'Reading file from {stats_root}')
            nat_fname = '{}/native_{}.de'.format(stats_root, feature)
            ref_fname = '{}/reference_{}.de'.format(stats_root, feature)
            if os.path.exists(nat_fname):
                nat_des += [DensityEstimate.read(nat_fname)]
            if os.path.exists(ref_fname):
                ref_des += [DensityEstimate.read(ref_fname)]

        nat_fname = '{}/native_{}.de'.format(merged, feature)
        ref_fname = '{}/reference_{}.de'.format(merged, feature)
        print(f'Write to {nat_fname}')
        DensityEstimate.merge(nat_des).write(nat_fname)
        print(f'Write to {ref_fname}')
        DensityEstimate.merge(ref_des).write(ref_fname)


# python -m src.combind.refit_stats merge-all /scratch/users/psuriana/combind_paper_systems/merged_stats/e3nn-run-20220422_101444-zybjdx87-best_val-max_helper_20 /scratch/users/psuriana/combind_paper_systems/stats/e3nn-run-20220422_101444-zybjdx87-best_val-max_helper_20

# TRAIN WITH TRAIN AND VAL SET
# python -m src.combind.refit_stats merge-all /scratch/users/psuriana/combind_paper_systems/merged_stats/e3nn-run-20220509_115644-29cuk25e-max_helper_20 /scratch/users/psuriana/combind_paper_systems/stats/e3nn-run-20220509_115644-29cuk25e-max_helper_20

# TRAIN WITH ALL TRAIN SET
# python -m src.combind.refit_stats merge-all /scratch/users/psuriana/combind_paper_systems/merged_stats/e3nn-run-20220508_160135-2kxspygq-max_helper_20 /scratch/users/psuriana/combind_paper_systems/stats/e3nn-run-20220508_160135-2kxspygq-max_helper_20

# GLIDE ONLY
# python -m src.combind.refit_stats merge-all /scratch/users/psuriana/combind_paper_systems/merged_stats/sp_es4-max_poses_300 /scratch/users/psuriana/combind_paper_systems/stats/sp_es4-max_poses_300
#
# python -m src.combind.refit_stats merge-all /scratch/users/psuriana/combind_paper_systems/merged_stats/e3nn_less_strict-run-20220422_171334-99eyoyz8-max_helper_20 /scratch/users/psuriana/combind_paper_systems/stats/e3nn_less_strict-run-20220422_171334-99eyoyz8-max_helper_20
@main.command()
@click.argument('merged_root')
@click.argument('protein_stats_root')
@click.option('--feature-names', default='hbond,saltbridge,contact,shape,mcss')
def merge_all(merged_root, protein_stats_root, feature_names):
    all_protein_dirs = set(glob(os.path.join(protein_stats_root, '*')))
    for protein_dir in all_protein_dirs:
        protein = Path(protein_dir).name
        merged_dir = os.path.join(merged_root, protein)
        os.makedirs(merged_dir, exist_ok=True)
        print(f'Write merged stats for {protein} to {merged_dir}')
        # Merge all other proteins except itself
        others = all_protein_dirs - set(protein_dir)
        merge_stats(merged_dir, list(others), feature_names)


if __name__ == '__main__':
    main()
