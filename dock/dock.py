from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.structutils.rmsd import ConformerRmsd
import os
import subprocess

GLIDE_ES4 = '''GRIDFILE  {grid}
LIGANDFILE   {ligands}
DOCKING_METHOD   confgen
POSES_PER_LIG   100
POSTDOCK_NPOSE   100
PRECISION   SP
NENHANCED_SAMPLING   4
'''

GLIDE_XP = '''GRIDFILE  {grid}
LIGANDFILE   {ligands}
DOCKING_METHOD   confgen
POSES_PER_LIG   100
POSTDOCK_NPOSE   100
PRECISION   XP
'''

GLIDE = '''GRIDFILE  {grid}
LIGANDFILE   {ligands}
DOCKING_METHOD   confgen
POSES_PER_LIG   30
POSTDOCK_NPOSE   30
PRECISION   SP
'''

GLIDE_ES4_REF = '''GRIDFILE  {grid}
LIGANDFILE   {ligands}
DOCKING_METHOD   confgen
POSES_PER_LIG   100
POSTDOCK_NPOSE   100
PRECISION   SP
NENHANCED_SAMPLING   4
USE_REF_LIGAND True
REF_LIGAND_FILE {reference}
CORE_DEFINITION mcssmarts
CORE_RESTRAIN   True
CORE_POS_MAX_RMSD 2
'''

GLIDE_REDOCK = '''GRIDFILE  {grid}
LIGANDFILE   {ligands}
DOCKING_METHOD   inplace
POSES_PER_LIG   100
POSTDOCK_NPOSE   100
PRECISION   SP
NENHANCED_SAMPLING   4
POSE_RMSD 1.
WRITE_RES_INTERACTION True
'''

def docking_failed(glide_log):
    if not os.path.exists(glide_log):
        return False
    with open(glide_log) as fp:
        logtxt = fp.read()
    phrases = ['** NO ACCEPTABLE LIGAND POSES WERE FOUND **',
               'NO VALID POSES AFTER MINIMIZATION: SKIPPING.',
               'No Ligand Poses were written to external file',
               'GLIDE WARNING: Skipping refinement, etc. because rough-score step failed.']
    return any(phrase in logtxt for phrase in phrases)

def dock(grid, ligands, root, name, enhanced,  reference, rescore, infile=None):
    if infile is None:
        infile = GLIDE_ES4 if enhanced else GLIDE
    print ("here "+ str(reference is None))
    print(rescore)
    if (not reference is None):
        infile = GLIDE_ES4_REF
    glide_in = '{}/{}.in'.format(root, name)
    glide_rescore_in = '{}/rescore-{}.in'.format(root, name)
    glide_pv = '{}/{}_pv.maegz'.format(root, name)
    glide_log = '{}/{}.log'.format(root, name)
    glide_cmd = 'glide -WAIT -LOCAL -RESTART {}'.format(os.path.basename(glide_in))

    if not rescore and os.path.exists(glide_pv):
        return

    if not rescore and enhanced and docking_failed(glide_log):
        return

    if not os.path.exists(root):
        os.system('mkdir {}'.format(root))
    if not rescore and not reference is None:

        with open(glide_in, 'w') as fp:
            fp.write(infile.format(grid=grid, ligands=ligands, reference=reference))
    elif rescore:
        print('dockign? ')
        with open(glide_rescore_in, 'w') as fp:
            fp.write(GLIDE_REDOCK.format(grid=grid, ligands=glide_pv))
        glide_cmd = 'glide -WAIT -LOCAL -RESTART {}'.format(os.path.basename(glide_rescore_in))

    else:
        with open(glide_in, 'w') as fp:
            fp.write(infile.format(grid=grid, ligands=ligands))
    print('help')
    subprocess.run(glide_cmd, cwd=root, shell=True)


#update
def atom_distance_and_features(native, pv):
    with StructureReader(native) as sts:
        native = list(sts)
        assert len(native) == 1, len(native)
        native = native[0]

    rmsds = []
    with StructureReader(pv) as reader:
        receptor = next(reader)
        for st in reader:
            conf_rmsd = ConformerRmsd(native, st)
            rmsds += [conf_rmsd.calculate()]
    return rmsds

def rmsd(native, pv):
    with StructureReader(native) as sts:
        native = list(sts)
        assert len(native) == 1, len(native)
        native = native[0]

    rmsds = []
    with StructureReader(pv) as reader:
        receptor = next(reader)
        for st in reader:
            conf_rmsd = ConformerRmsd(native, st)
            rmsds += [conf_rmsd.calculate()]
    return rmsds

def filter_native(native, pv, out, thresh):
    with StructureReader(native) as sts:
        native = list(sts)
        assert len(native) == 1, len(native)
        native = native[0]

    near_native = []
    with StructureReader(pv) as reader:
        receptor = next(reader)
        for st in reader:
            conf_rmsd = ConformerRmsd(native, st)
            if conf_rmsd.calculate() < thresh:
                near_native += [st]

    print('Found {} near-native poses'.format(len(near_native)))
    if not near_native:
        print('Resorting to native pose.')
        native.property['r_i_docking_score'] = -10.0
        near_native = [native]

    with StructureWriter(out) as writer:
        writer.append(receptor)
        for st in near_native:
            writer.append(st)
