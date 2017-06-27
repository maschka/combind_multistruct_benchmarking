#!/share/PI/rondror/software/miniconda/bin/python
import os
import sys
import time

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1/run"
FUZZY_SCRIPT = "/share/PI/rondror/jetynan/combind/2_fingerprint/fuzzyifp.py"

DATA = "/scratch/PI/rondror/docking_data/"
dataset = sys.argv[1]

os.chdir(DATA + dataset)

output_dir = 'structure_fingerprints2'
OUTPUT = DATA + dataset + '/' + output_dir + '/ifp.fp'

processed = [o for o in os.listdir("./processed/") if os.path.isfile(os.path.join("./processed/",o))]
ligands = [o for o in os.listdir("./ligands/") if os.path.isfile(os.path.join("./ligands/",o))]

if not os.path.exists(output_dir):
    os.system("mkdir " + output_dir)
os.chdir(output_dir)

for structFile in processed:
    struct = os.path.splitext(structFile)[0]
    ligandFile = filter(lambda x: struct.lower() in x.lower(), ligands)[0]
    with open(struct + '.sh', 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('OUTPUT=$(' + SCHRODINGER + ' ' + FUZZY_SCRIPT + ' -receptor ' + DATA + dataset + "/processed/" + structFile + ' -ligand ' + DATA + dataset + "/ligands/" + ligandFile + ' -input pair)\n')
        f.write("echo \"" + struct + ";" + '$OUTPUT\" >> ' + OUTPUT)
        
    os.system("sbatch --time=10:00:00 -n 1 -p rondror " + struct + '.sh')
