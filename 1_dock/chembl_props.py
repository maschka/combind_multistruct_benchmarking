import os
import sys

def read_molw(dir_path=None):
    molw = {}
    fpath = 'chembl/molw.txt'
    if dir_path is not None:
        fpath = '{}/{}'.format(dir_path, fpath)
    if not os.path.exists(fpath):# and os.path.exists('chembl'):
        print('missing molw file')
        return {}
    with open(fpath) as fname:
        for line in fname:
            chembl, mw = line.strip().split(',')
            molw[chembl] = float(mw)
    return molw

def write_props(lm):
    if len(lm.all_ligs) == 0: return
    from schrodinger.structure import StructureReader
    l_path = 'ligands/prepared_ligands/{}/{}.mae'
    write_molw(lm,l_path,StructureReader)
    write_duplicates(lm,l_path,StructureReader)

def write_molw(lm,l_path,StructureReader):
    with open('chembl/molw.txt', 'w') as fname:
        for f in lm.all_ligs:
            mw = next(StructureReader(l_path.format(f,f))).total_weight
            fname.write('{},{}\n'.format(f,mw))

def write_duplicates(lm,l_path,StructureReader):
    duplicates = {}
    all_ligs = {l : next(StructureReader(l_path.format(l,l))) for l in lm.all_ligs}
    for l1 in all_ligs:
        for l2 in duplicates:
            if all_ligs[l1].isEquivalent(all_ligs[l2], True):
                duplicates[l2].append(l1)
                break
        else:
            duplicates[l1] = [l1]
    with open('chembl/duplicates.txt','w') as f:
        for lig,l_list in duplicates.items():
            f.write('{}\n'.format(','.join(l_list)))

def read_duplicates(dir_path=None):
    duplicates = {}
    unique = set([])
    fpath = 'chembl/duplicates.txt'
    if dir_path is not None:
        fpath = '{}/{}'.format(dir_path, fpath)
    if not os.path.exists(fpath): 
        print('missing duplicates file')
        return unique, duplicates
    with open(fpath) as f:
        for line in f:
            l_list = line.strip().split(',')
            if len(l_list) == 1: unique.add(l_list[0])
            else:
                duplicates[l_list[0]] = set(l_list)
    return unique, duplicates







