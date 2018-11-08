import sys
import os
import numpy as np
from density_estimate import DensityEstimate
from pairs import LigPair
from containers import Protein
from shared_paths import shared_paths, feature_defs

def merge_stats(stats1, stats2, weight):
    '''
    stats {'native' or 'reference': {interaction (str): DensityEstimate}}
    
    Merges the dictionaries of DensityEstimate's stats1 and stats2.
    '''
    out = {}
    distributions = set(list(stats1.keys()) + list(stats2.keys()))
    for d in distributions:
        out[d] = {}
        interactions = []
        if d in stats1: interactions += list(stats1[d].keys())
        if d in stats2: interactions += list(stats2[d].keys())
        for i in set(interactions):
            if   d not in stats1 or i not in stats1[d]:
                out[d][i] = stats2[d][i]
            elif d not in stats2 or i not in stats2[d]:
                out[d][i] = stats1[d][i]
            else:
                out[d][i] = stats1[d][i].average(stats2[d][i], weight)
    return out

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

def weighting(w_ref, protein):
    '''
    w_ref [(float, float), ...]: pairs of glide scores
    protein (str): name of protein

    Returns glide scores mapped to probabilities using the method specified
    in shared_paths['stats']['weighting'].
    '''
    pnative = '{}/{}/stats/{}/pnative.de'.format(shared_paths['data'],
                                            protein,
                                            shared_paths['stats']['version'])
    if shared_paths['stats']['weighting'] == 'absolute':
        pnative = DensityEstimate.read(pnative)
        return np.array([pnative(g1) * pnative(g2)
                         for (g1, g2) in w_ref])
    elif shared_paths['stats']['weighting'] == 'relative':
        pnative = DensityEstimate.read(pnative)
        g1_top, g2_top = lig_pair.get_gscores(0, 0)
        return np.array([pnative(g1-g1_top) * pnative(g2-g2_top)
                         for (g1, g2) in w_ref])
    elif shared_paths['stats']['weighting'] == 'unweighted':
        return 1
    else:
        assert False

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
            w_ref += [lig_pair.get_gscores(r1, r2)]

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
        w_ref = weighting(w_ref, protein)

        for d, X, w in [('native', X_native, 1), ('reference', X_ref, w_ref)]:
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
    lm = Protein(protein).lm
    ligands = lm.docked(lm.pdb)[:shared_paths['stats']['n_ligs']+1]
    self_docked = lm.st+'_lig'
    if self_docked in ligands:
        ligands.remove(self_docked)
    else:
        ligands.pop(-1)

    protein_stats = {}
    for j, ligand1 in enumerate(ligands):
        for ligand2 in ligands[j+1:]:
            ligand_stats = statistics_lig_pair(protein, ligand1, ligand2,
                                               interactions)
            stats = merge_stats(stats, ligand_stats, shared_paths['stats']['ligands_equal'])
    for d, interactions in stats.items():
        for i, de in interactions.items():
            de.write(fname.format(i, d))
    return protein_stats

def statistics(data, interactions):
    '''
    data    {protein (str): [ligand (str), ...]}
    interactions {name (str): [code (int), ...]}

    Returns DensityEstimate's representing the native and reference
    distribution of the statistics for all proteins, ligands in data.
    '''
    stats = {}
    for protein in data:
        protein_stats = statistics_protein(protein, interactions)
        stats = merge_stats(stats, protein_stats, shared_paths['stats']['proteins_equal'])
    return stats

def gscore_statistics(data):
    '''
    data    {protein (str): [ligand (str), ...]}

    Returns DensityEstimate's representing the native and reference
    distribution of the glide scores for all proteins, ligands in data.
    '''
    stats = {}
    params = shared_paths['stats']
    for protein, ligands in data.items():
        protein_stats = {}
        prot = Protein(protein)
        prot.load_docking(ligands, load_fp=False, load_mcss=False)
        st = prot.lm.st
        docking = prot.docking[st]
        for ligand in ligands:
            # Get input data
            poses = docking.ligands[ligand].poses
            native = np.array([pose.rmsd <= params['native_thresh']
                               for pose in poses[:params['max_poses']]])
            glide  = np.array([pose.gscore
                               for pose in poses[:params['max_poses']]])

            if glide.shape[0] and params['weighting'] == 'relative':
                glide -= glide.min()
            
            # Compute densities
            ligand_stats = {}
            for k, X in [('native', glide[native==1]), ('reference', glide)]:
                ligand_stats[k] = {}
                ligand_stats[k][0] = DensityEstimate(
                    points  = params['gscore_points'],
                    reflect = params['weighting'] == 'relative',
                    sd      = params['gscore_sd'],
                    domain  = params['gscore_domain'],
                    out_of_bounds = 0)
                ligand_stats[k][0].fit(X)
            
            protein_stats = merge_stats(protein_stats, ligand_stats)
        stats = merge_stats(stats, protein_stats)
    return stats['native'][0], stats['reference'][0]


if __name__ == '__main__':
    '''
    python statistics.py mode protein

    If mode is 'gscore', compute probability native given gscore for all
    proteins other than 'protein'.
    
    Else, compute overlap statistics for ligand1 and ligand2 for interactions,
    which should be passed as a comma separated list.
    '''
    mode, protein = sys.argv[1:]
    if mode == 'gscore':
        protein = sys.argv[2]
        os.chdir(shared_paths['data'])
        datasets = [d for d in sorted(os.listdir('.'))
                    if d[0] != '.' and d != protein]
        data = {}
        for d in datasets:
            lm = Protein(d).lm
            data[d] = lm.docked(lm.pdb)[:shared_paths['stats']['n_ligs']+1]
            self_docked = lm.st+'_lig'
            if self_docked in data[d]:
                data[d].remove(self_docked)
            else:
                data[d].pop(-1)
        native, reference = gscore_statistics(data)
        fname = '{}/{}/stats/{}/{}.de'.format(shared_paths['data'],
                                              protein,
                                              shared_paths['stats']['version'],
                                              '{}')
        native.write(fname.format('native'))
        reference.write(fname.format('reference'))
        native.ratio(reference, prob = False).write(fname.format('pnative'))

    elif mode == 'fp':  
        statistics([protein], feature_defs.keys())
    else:
        assert False
