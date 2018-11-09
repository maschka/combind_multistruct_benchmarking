import numpy as np
from density_estimate import DensityEstimate
from pairs import LigPair
from containers import Protein
from shared_paths import shared_paths, feature_defs

# Helper Functions
def merge_dicts_of_lists(d1, d2):
    assert (type(d1) == dict) == (type(d2) == dict)
    if type(d1) == dict:
        return {key: merge_dicts_of_lists(d1[key] if key in d1 else [],
                                          d2[key] if key in d2 else [])
                for key in set(d1.keys()).union(set(d2.keys()))}
    if type(d1) != list: d1 = [d1]
    if type(d2) != list: d2 = [d2]
    return d1 + d2

def merge_stats(stats, weight):
    for d, interactions in stats.items():
        for i, des in interactions.items():
            stats[d][i] = DensityEstimate.merge(des, weight)
    return stats

# I/O
def get_fname(protein, ligand1, ligand2):
    '''
    protein, ligand1, ligand2 (str): names of protein and ligands
    for protein level stats, set ligand1 and ligand2 to None.

    Returns path to where statistics files should be read/written.
    '''
    assert (ligand1 is None) == (ligand2 is None)
    if ligand1 is None and ligand2 is None:
        ID = protein
    else:
        ID = '{}-{}'.format(ligand1, ligand2)
    version = shared_paths['stats']['version']
    return "{}/{}/stats/{}/{}-{}-{}.de".format(shared_paths['data'],
                                               protein, version, ID,
                                               '{}', '{}')

def read_stats(fname, interactions):
    '''
    fname (str): output from get_fname(...)
    interactions [interaction (str),  ...]: interactions to read from files.
    '''
    stats = {'native':{}, 'reference':{}}
    for d in stats:
        for i in interactions:
            try: 
                stats[d][i] = DensityEstimate.read(fname.format(i, d))
            except:
                pass
    return stats

# Core stats computation
def get_interaction(lig_pair, interaction):
    '''
    lig_pair (LigPair): ligand pair for which to get features
    interaction (str): interaction to get scores for
    
    Returns X_native: features for native poses
            X_ref:    features for all poses
            w_ref:    pairs of glide scores for all poses
    '''
    X_native, X_ref, w_ref = [], [], []
    for (r1,r2), pp in lig_pair.pose_pairs.items():
        if max(r1, r2) > shared_paths['stats']['max_poses']: continue
        pp_x = lig_pair.get_feature(interaction, r1, r2)
        if pp_x is not None:
            if pp.correct():
                X_native += [pp_x]
            X_ref += [pp_x]

    X_native, X_ref = np.array(X_native), np.array(X_ref)
    return X_native, X_ref,  w_ref

def statistics_lig_pair(protein, ligand1, ligand2, interactions):
    '''
    protein (str):  protein name
    ligand1,2 (str): names of ligands for which to compute statistics
    interactions [interaction (str), ...]: interactions for which to compute
        statistics
    pnative (function): Maps glide scores to probability correct. Use a
        DensityEstimate if you wish to weight by glide score, else pass
        lambda x: 1
    sd (float): standard deviation for density estimate.
    '''
    # If all statistics have already been computed, just read and return.
    fname = get_fname(protein, ligand1, ligand2)
    stats = read_stats(fname, interactions)
    if all(    interaction in stats['native']
           and interaction in stats['reference']
           for interaction in interactions):
        return stats
    
    print('Computing statistics for:', protein, ligand1, ligand2)
    # Load relevant data.
    prot = Protein(protein)
    prot.load_docking([ligand1, ligand2],
                      load_fp = True, load_mcss = 'mcss' in interactions)
    lm = prot.lm
    docking = prot.docking[lm.st]
    lig_pair = LigPair(docking.ligands[ligand1],
                       docking.ligands[ligand2],
                       interactions, lm.mcss if 'mcss' in interactions else None,
                       shared_paths['stats']['max_poses'])

    # Compute all remaining statistics and write to files.
    for interaction in interactions:
        if (    interaction in stats['native']
            and interaction in stats['reference']):
            continue

        X_native, X_ref, w_ref = get_interaction(lig_pair, interaction)
        w_ref = 1

        for d, X in [('native', X_native), ('reference', X_ref)]:
            domain = (0, 15) if interaction == 'mcss' else (0, 1)

            stats[d][interaction] = DensityEstimate(domain = domain,
                              sd = shared_paths['stats']['stats_sd']*(domain[1]-domain[0]),
                              reflect = True)
            stats[d][interaction].fit(X)
            stats[d][interaction].write(fname.format(interaction, d))
    return stats

def statistics_protein(protein, interactions):
    fname = get_fname(protein, None, None)
    stats = read_stats(fname, interactions)
    if all(    i in stats['native']
           and i in stats['reference']
           for i in interactions):
        return stats
    
    print('Computing statistics for:', protein)
    ligands = Protein(protein).get_xdocked_ligands(shared_paths['stats']['n_ligs'])
    print(ligands)

    for j, ligand1 in enumerate(ligands):
        for ligand2 in ligands[j+1:]:
            ligand_stats = statistics_lig_pair(protein, ligand1, ligand2, interactions)
            stats = merge_dicts_of_lists(stats, ligand_stats)

    stats = merge_stats(stats, shared_paths['stats']['poses_equal'])
    
    for d, interactions in stats.items():
        for i, de in interactions.items():
            de.write(fname.format(i, d))
    return stats

def statistics(data, interactions):
    '''
    data    {protein (str): [ligand (str), ...]}
    interactions {name (str): [code (int), ...]}

    Returns DensityEstimate's representing the native and reference
    distribution of the statistics for all proteins, ligands in data.
    '''
    stats = {'native': {}, 'reference':{}}
    for protein in data:
        protein_stats = statistics_protein(protein, interactions)
        stats = merge_dicts_of_lists(stats, protein_stats)
    return merge_stats(stats, shared_paths['stats']['ligands_equal'])

if __name__ == '__main__':
    import sys
    protein = sys.argv[1]
    statistics([protein], feature_defs.keys())
