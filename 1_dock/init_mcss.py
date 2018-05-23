import os
import sys
import itertools

from parse_chembl import load_chembl_proc

from schrodinger.structure import StructureReader, StructureWriter, _StructureProperty

atom_typing = 7
out_dir = 'ligands/mcss/mcss{}'.format(atom_typing)

command = '$SCHRODINGER/utilities/canvasMCS -JOB {} -WAIT -imae {}_in.mae -omae {}.mae -atomtype {} -nobreakaring -exclusive'

l_dir = 'ligands/prepared_ligands'
l_path = '{}/{}/{}.mae'

group_size=100
def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n 
    return itertools.izip_longest(*args, fillvalue=fillvalue)

def init_mcss(chembl=None, max_num=25):

    all_pairs = set([])
    if chembl is not None:
        for f, f_data in chembl.items():
            for q,c in f_data.items():
                for i in range(len(c)):
                    for j in range(i+1,len(c)):
                        all_pairs.add((c[i],c[j]))
    else:
        chembl_info = load_chembl_proc()
        all_ligs = sorted([l for l in os.listdir(l_dir) if os.path.exists(l_path.format(l_dir,l,l))])
        pdb_ligs = [l for l in all_ligs if l[:6] != 'CHEMBL']
        chembl_ligs = [l for l in all_ligs if l in chembl_info and chembl_info[l].ki <= 1000 and chembl_info[l].mw <= 1000]
        chembl_ligs.sort(key=lambda x: chembl_info[x].ki)

        os.system('mkdir -p {}'.format(out_dir))

        for i1 in range(min(max_num,len(pdb_ligs))):
            for i2 in range(i1+1,min(max_num,len(pdb_ligs))):
                l1, l2 = pdb_ligs[i1], pdb_ligs[i2]
                all_pairs.add((l1, l2))

        for l1 in pdb_ligs:
            for l2 in chembl_ligs:
                all_pairs.add((l1, l2))

    mcss_compute(all_pairs)

def mcss_compute(all_pairs):
    no_mcss = []
    no_size_file = []

    for l1,l2 in all_pairs:
        if not os.path.exists(l_path.format(l_dir,l1,l1)): continue
        if not os.path.exists(l_path.format(l_dir,l2,l2)): continue
        if os.path.exists('{}/{}-{}.txt'.format(out_dir, l1, l2)):
            if os.path.exists('{}/{}-{}.mae'.format(out_dir, l1, l2)):
                # delete the structure file once the text file is written
                os.system('rm -f {}/{}-{}.mae'.format(out_dir, l1, l2))
                os.system('rm -f {}/{}-{}_in.mae'.format(out_dir, l1, l2))

        elif not os.path.exists('{}/{}-{}.mae'.format(out_dir, l1, l2)):
            no_mcss.append((l1,l2))
        elif not os.path.exists('{}/{}-{}.txt'.format(out_dir, l1, l2)):
            no_size_file.append((l1,l2))

    if len(no_mcss) > 0:
        print len(no_mcss), 'mcss pairs left'
        get_mcss(no_mcss)
    if len(no_size_file) > 0:
        print len(no_size_file), 'size pairs left'
        get_size_file(no_size_file)

def get_mcss(mcss_pairs):
    for i, pairs in enumerate(grouper(group_size, mcss_pairs)):
        with open('{}/mcss{}.sh'.format(out_dir,i), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            for p in pairs:
                if p is None: continue
                l1,l2 = p
                name = '{}-{}'.format(l1, l2) 
                st1 = StructureReader(l_path.format(l_dir,l1,l1)).next()
                st2 = StructureReader(l_path.format(l_dir,l2,l2)).next()
                st1.title = l1
                st2.title = l2

                st_wr = StructureWriter('{}/{}_in.mae'.format(out_dir,name))
                st_wr.append(st1)
                st_wr.append(st2)
                st_wr.close()
                    
                f.write(command.format(name, name, name, atom_typing)+'\n')
            f.write('wait\n')
        os.chdir(out_dir)
        os.system('sbatch -p owners --tasks=6 --cpus-per-task=1 -t 0:30:00 -o mcss.out mcss{}.sh'.format(i))
        os.chdir('../../..')

def get_size_file(size_pairs):
    for l1,l2 in size_pairs:
        with open('{}/{}-{}.txt'.format(out_dir, l1, l2),'w') as f:
            for st in StructureReader('{}/{}-{}.mae'.format(out_dir, l1, l2)):
                prop = _StructureProperty(st)
                smarts = prop['s_canvas_MCS_SMARTS']
                msize = prop['i_canvas_MCS_Atom_Count'] 
                size = len([a for a in st.atom if a.element != 'H'])

                f.write('{},{},{},{}\n'.format(st.title, size, msize, smarts))



