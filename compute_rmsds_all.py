#!/bin/env python

import pandas as pd
import numpy as np
import click
import os
from glob import glob
from schrodinger.structure import StructureReader, StructureWriter
from utils import *

###############################################################################

# Defaults
stats_root = os.environ['COMBINDHOME'] + '/stats_data/default'
mcss_version = 'mcss16'
shape_version = 'pharm_max'
ifp_version = 'rd1'


# @click.group()
# def main():
#     pass


from subprocess import run

command = ('$SCHRODINGER/utilities/structalign '
           '-asl "(not chain. L and not atom.element H) and (fillres within {0} chain. L)" '
           '-asl_mobile "(not chain. L and not atom.element H) and (fillres within {0} chain. L)" '
           '{1} {2}')

import os
from schrodinger.structure import StructureReader


def split_complex(st, pdb_id, template):
    os.system('mkdir -p structures/proteins structures/ligands')
    lig_path = 'structures/ligands/{}_lig_to_{}.mae'.format(pdb_id, template)
    prot_path = 'structures/proteins/{}_prot_to_{}.mae'.format(pdb_id, template)

    if not os.path.exists(lig_path) and len([a.index for a in st.atom if a.chain == 'L']) > 0:
        lig_st = st.extract([a.index for a in st.atom if a.chain == 'L'])
        lig_st.title = '{}_lig_to_{}'.format(pdb_id, template)
        lig_st.write(lig_path)

    if not os.path.exists(prot_path):
        prot_st = st.extract([a.index for a in st.atom if a.chain != 'L'])
        prot_st.title = '{}_prot_to_{}'.format(pdb_id, template)
        prot_st.write(prot_path)


def struct_sort(structs):
    for struct in structs:
        for template in structs:
            opt_complex = 'structures/aligned/{}/rot-{}_query_to_{}.mae'.format(struct, struct, template)

            if os.path.exists(opt_complex):
                comp_st = next(StructureReader(opt_complex))
                split_complex(comp_st, struct, template)


def align_successful_all(out_dir, struct, template):
    if not os.path.exists('{}/{}/rot-{}_query_to_{}.mae'.format(out_dir, struct, struct, template)):
        return False

    # if os.path.exists('{}/{}/{}_template.mae'.format(out_dir, struct, struct)):
    # return True # query = template so we don't need to check alignment

    with open('{}/{}/align_to_{}.out'.format(out_dir, struct, template), 'r') as f:
        for line in f:
            tmp = line.strip().split()
            if len(tmp) > 0 and tmp[0] == 'Alignment':
                if float(tmp[2]) > 0.4:
                    print('-- Alignment warning!', struct, float(tmp[2]))
                    return False
                return True
        else:
            print('alignment failure', struct)
            return False


# @main.command()
# @click.argument('docking', default='')
# @click.argument('crystal', default='')
# changed crystal path - 11/14 - john
def rmsd_all(docking='docking/*/*_pv.maegz', crystal='structures/ligands/*_lig_to_*.mae'):
    """
    Compute rmsd of docked poses to a reference pose.

    docking and crystal should be paths to the docked poses in pose viewer
    format and the reference poses. Non-expanded patterns capturing multiple
    docking results and references (for different ligands) can be provided to
    process multiple files at a time.

    It is required that the docked poses and crystal ligands have the same
    name. For the docking the name is the part of the basename before '-to-' and
    for the reference poses it is the part of the basename before '_lig'.
    """
    from dock.dock import rmsd
    print(crystal)
    docking = glob(docking)
    crystal = glob(crystal)

    def docking_to_name(path):
        path = path.split('/')[-1]
        ligand = path.split('-to-')[0]
        ligand = ligand.split('_lig')[0]
        struct = path.split('-to-')[1].split('_pv.')[0]
        return ligand, struct

    def crystal_to_name(path):
        path = path.split('/')[-1]
        ligand = path.split('_lig_to_')[0]
        template = path.split('_lig_to_')[1].split('.')[0]
        return ligand, template

    for docking_path in docking:
        # if os.path.exists(out):
        # continue

        name = docking_to_name(docking_path)
        print(name)
        # print()
        crystal_path = [crystal_path for crystal_path in crystal
                        if crystal_to_name(crystal_path) == name]
        print(crystal_path)
        if len(crystal_path) == 0:
            print('No crystal pose for {}: {}'.format(name, docking_path))
            continue
        if len(crystal_path) > 1:
            print('Multiple crystal poses for {}: {}. Doing nothing.')
            continue

        crystal_path = crystal_path[0]
        print('Computing rmsd for {} to {}.'.format(docking_path, crystal_path))
        for path in crystal_path:
            out = docking_path.replace('.maegz', '_rmsd_to.npy')
            rmsds = rmsd(crystal_path, docking_path)
        np.save(out, rmsds)


def struct_align_all(template, structs, dist=15.0, retry=True,
                     # processed_out='structures/processed/{pdb}/{pdb}_out.mae',
                     processed_out='structures/processed/{pdb}/{pdb}_merged.mae',
                     align_dir='structures/aligned'):
    # template_path = 'structures/aligned/{}/{}_query.mae'.format(template, template)
    # tempalte_path = 'structures/processed/{}/{}_merged.mae'.format(template, template)
    template_path = processed_out.format(pdb=template)
    if not os.path.exists(template_path):
        print('template not processed', template_path)
        return
    for struct in structs:
        query_path = processed_out.format(pdb=struct)
        # if not os.path.exists(query_path) or
        if align_successful_all(align_dir, struct, template):
            continue

        print('align', struct, template)

        os.system('mkdir -p {}'.format(align_dir))
        # os.system('rm -rf {}/{}'.format(align_dir, struct))
        os.system('mkdir -p {}/{}'.format(align_dir, struct))

        _workdir = '{}/{}'.format(align_dir, struct)
        _template_fname = '{}_template.mae'.format(template)
        _query_fname = '{}_query_to_{}.mae'.format(struct, template)

        os.system('cp {} {}/{}'.format(template_path, _workdir, _template_fname))
        os.system('cp {} {}/{}'.format(query_path, _workdir, _query_fname))

        with open('{}/align_in_to_{}.sh'.format(_workdir, template), 'w') as f:
            f.write(command.format(dist, _template_fname, _query_fname))
        run('sh align_in_to_{}.sh > align_to_{}.out'.format(template, template), shell=True, cwd=_workdir)

        if retry and not align_successful_all(align_dir, struct, template):
            print('Alignment failed. Trying again with a larger radius.')
            struct_align_all(template, [struct], dist=25.0, retry=False,
                             processed_out=processed_out, align_dir=align_dir)
