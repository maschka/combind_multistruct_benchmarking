import os
import sys

from shared_paths import shared_paths
from sort_downloads import sort_downloads

from align_structs import align_structs
from process_structs import process_structs
from sort_files import sort_files

from grids import make_grids
from dock import dock
from fp import fp
from MCSSController import compute_mcss, verify_mcss
from stats import stats

from chembl_sort import get_ligands, proc_ligands
from chembl_props import write_props
from pick_helpers import pick_helpers, load_helpers
from score import score

from verify_docking import check_docked_ligands

sys.path.append('../3_analyze')
from containers import LigandManager

os.chdir(shared_paths['data'])

todo = list(sys.argv[1])
if len(todo) == 0:
    todo = list('12345') 

datasets = sys.argv[2:]
if datasets == []:
    datasets = [d for d in sorted(os.listdir('.')) if d[0] != '.' and d[-3:] != 'old']

datasets=reversed(datasets)
for i, d in enumerate(datasets):
    print(d, i)

    os.chdir(d)
    
    lm = LigandManager(shared_paths, d)

    if '0' in todo:
        sort_downloads()
 
    # 1. prepare proteins and associated grids 
    if '1' in todo:
        process_structs()             # Runs prepwizard
        align_structs()               # Align and give consistent numbering
        sort_files()                  # Creates ligand, protein, and complex directories
        make_grids()                  # Creates grid for all proteins
        dock(lm)                      # Dock PDB ligands to lm.st
        dock(lm, mode = 'inplace')
        dock(lm, mode = 'mininplace')
        dock(lm, mode = 'XP')
        dock(lm, mode = 'expanded')
     
    # 2. prepare ligands
    if '2' in todo:
       get_ligands()                  # Writes MAE files for all ligs to ligands/raw_files
       proc_ligands()                 # Runs prepwizard & epik on all ligs
       
    if 'm' in todo:
        compute_mcss(lm, compute_rmsds = False) # Computes MCSS, for use in pick_helpers

    if 'v' in todo:
        verify_mcss(lm)
    if 'g' in todo:
        check_docked_ligands(lm)
    # force redo of chembl info (do this if new chembl ligands have been added)
    # Do this after all MCSS files have been written!
    if 'c' in todo: #pass
         os.system('rm chembl/helpers/*')
         os.system('rm chembl/duplicates.txt')
         os.system('rm chembl/molw.txt')
         os.system('rm chembl/macrocycle.txt') 
         write_props(lm)

    # # 3. decide what ligands to use and prepare them
    if '3' in todo:
        pick_helpers(lm)              # Picks chembl ligands for use in scoring for each pdb ligand
        dock(lm, load_helpers())      # Dock chembl ligands to be used for scoring all pdb ligands
        fp(lm)                        # Writeout fingerprints for docked poses and pdb structures
        compute_mcss(lm, load_helpers())  # Performs all phases of MCSS computation

    # 4. compute statistics
    if '4' in todo:
        stats(lm)                     # Computes statistics and writes to stats/stats11

    # 5. run method on all ligands
    if '5' in todo:
        score(lm, load_helpers())     # Computes optimal cluster

    os.chdir('..')
