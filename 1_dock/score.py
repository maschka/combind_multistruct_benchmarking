import os
import sys
import itertools

group_size=10
def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n 
    return itertools.izip_longest(*args, fillvalue=fillvalue)

output_dir = 'scores/scores10'
cmd = '$SCHRODINGER/run /scratch/PI/rondror/jbelk/method/combind/3_analyze/scores.py {} {} {}'

settings = {
    'data_dir' : '/scratch/PI/rondror/jbelk/method/data',
    'glide_dir' : 'docking/glide12',
    'ifp_dir' : 'ifp/ifp3',
    'mcss_dir' : 'mcss/mcss7',
    'stats_dir': 'stats/stats2',
    'k_list' : ['mcss','hbond','sb1'],#,'contact'],#,'pipi','contact']
    'num_stats_ligs' : 10,
    'normalize' : True,
    'num_pred_chembl' : 10,
    'num_poses' : 100,
    't' : 10,
    'mcss_sort': False, # set to True when using best_mcss.txt
    'chembl_file': 'best_affinity.txt',
    'score_mode': 'ALL'
    #'use_chembl':False
}

#stats_exclude = [ set(['B1AR','B2AR']), set(['AR','ERA']) ]

def write_settings_file(out_path, settings):
    #if os.path.exists(out_path): return
    with open(out_path,'w') as f:
        for varname, var in settings.items():
            if type(var) is str: var = '"{}"'.format(var)
            f.write('{}={}\n'.format(varname, var))

def score(lm, helpers):
    all_p = [d for d in sorted(os.listdir(settings['data_dir'])) if d[0] != '.' and d[-3:] != 'old']
    settings['stats_prots'] = [p for p in all_p if p != lm.prot]
    os.system('mkdir -p {}'.format(output_dir))
    os.chdir(output_dir)
    #os.system('rm -f *')
    write_settings_file('settings.py', settings)

    unfinished = sorted([l for l in helpers[settings['chembl_file']] 
        if not os.path.exists('{}-to-{}.sc'.format(l,lm.st))])

    if len(unfinished) > 0:
        print len(unfinished), 'scores left'

    for i,group in enumerate(grouper(group_size, unfinished)):
        with open('{}.sh'.format(i),'w') as f:
            f.write('#!/bin/bash\n')
            f.write(cmd.format(lm.st, lm.prot, ' '.join([q for q in group if q is not None])))
        os.system('sbatch -t 1:00:00 -p rondror {}.sh'.format(i))
        #break
    os.chdir('../..')



