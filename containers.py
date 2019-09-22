"""
These classes keep track of
Poses, Ligands, Docking results for a given grid,
ligand names and preparation progress (LigandManager),
and ligands + docking results for a Protein.

Design decisions:
 - Parameterize by passing directory paths
   - Exception for struct, since needs to be specified
     in the docking file name
 - Open all files with absolute paths
 - Only load data if explicitly specified
 - Have default structure specified by LigandManager.
   always use this structure unless explicitly specified
   otherwise.
"""

import os
import sys
import numpy as np

from shared_paths import shared_paths
from ifp.fp_controller import parse_fp_file
from dock.parse_chembl import load_chembl_proc
from dock.chembl_props import read_duplicates
from mcss.mcss_controller import MCSSController

class Pose:
    """
    Represents a single pose of a ligand.

    Parameterized by (protein, structure, ligand, rank).
    """
    def __init__(self, rmsd, gscore, emodel, fp):
        """
        rmsd (float): RMSD to crystallographic pose
        gscore (float): glide score
        emodel (float): emodel (from glide)
        fp ({(int, string): float}): dict mapping interactiontype, resname pairs to interaction score
        """
        self.rmsd = rmsd
        self.gscore = gscore
        self.emodel = emodel
        self.fp = fp

class Ligand:
    """
    Stores all poses of a given ligand docked to a given structure.

    This class defines file path convention for docking and ifp subdirectories.

    Parameterized by (protein, structure, ligand).
    """
    def __init__(self, ligand, dock_dir, fp_dir, struct):
        self.ligand = ligand
        self.glide_path = '{}/{}-to-{}'.format(dock_dir, ligand, struct)
        self.fp_path = '{}/{}-to-{}-{}.fp'.format(fp_dir, ligand, struct,
                                                  shared_paths['docking'])
        self.crystal_fp_path = '{}/{}.fp'.format(fp_dir,
                                                 ligand.replace('_crystal_lig', '_struct'))

        self.poses = None

    def load_poses(self, load_fp):
        gscores, emodels, rmsds = self.parse_glide_output()

        fps = {}
        if load_fp:  
            fps = parse_fp_file(self.fp_path)
            assert len(fps) >= min(shared_paths['stats']['max_poses'], len(rmsds)), \
                   ('missing for {}'.format(self.fp_path))
 
        assert len(gscores) == len(rmsds) == len(emodels), \
               ('Number of glide scores, rmsds, and emodels '
                'is not equal for {}'.format(self.glide_path))
        assert gscores == sorted(gscores), \
               'Glide scores are not correctly ordered for {}'.format(self.glide_path)

        self.poses = [Pose(rmsds[i], gscores[i], emodels[i], fps.get(i, {}))
                      for i in range(len(gscores))]

    def load_crystal_pose(self):
        try:
            fps = parse_fp_file(self.crystal_fp_path)
            self.poses = [Pose(0, -1000, -1000, fps[0])]
        except IOError:
            pass

    def parse_glide_output(self):
        if not os.path.exists(self.glide_path):
            return [], [], []
        pair = self.glide_path.split('/')[-1]
        if os.path.exists('{}/{}_pv.maegz'.format(self.glide_path, pair)):
            return self.parse_rept_file(pair)
        else:
            print('not finished', self.glide_path)
            return [], [], []

    def parse_rept_file(self, pair):
        lig, prot = pair.split('-to-')
        rept_file = '{}/{}.rept'.format(self.glide_path, pair)
        rmsd_file = '{}/rmsd.csv'.format(self.glide_path, pair)
        
        gscores, emodels, rmsds = [], [], []
        with open(rept_file) as fp:
            for line in fp:
                line = line.strip().split()
                if len(line) <= 1 or (line[1] != lig and line[1] != lig+'_out' and line[1] != '1'): continue
                rank, lig_name, lig_index, score = line[:4]
                emodel = line[13]
                if line[1] == '1':
                    rank, lig_index, score = line[:3]
                    emodel = line[12]
                gscores.append(float(score))
                emodels.append(float(emodel))
            
        if not os.path.exists(rmsd_file):
            return gscores, emodels, [None]*len(gscores)

        with open(rmsd_file) as fp:
            for line in fp:
                line = line.strip().split(',')
                if line[3] == '"RMSD"': continue
                rmsds.append(float(line[3][1:-1]))

        return gscores, emodels, rmsds

class Docking:
    """
    Stores docking results for a set of ligands to a single structure.

    This class defines path to docking and ifp directories.

    Parameterized by (protein, structure).
    """
    def __init__(self, root, struct):
        self.dock_dir = '{}/docking/{}'.format(root, shared_paths['docking'])
        self.ifp_dir = '{}/ifp/{}'.format(root, shared_paths['ifp']['version'])
        self.struct = struct
        self.ligands = {}
        self.num_poses = {}

    def load(self, ligands, load_fp, load_crystal):
        for ligand in ligands:
            pair = '{}-to-{}'.format(ligand.replace('_crystal', ''), self.struct)
            self.ligands[ligand] = Ligand(ligand, self.dock_dir, self.ifp_dir, self.struct)
            if load_crystal:
                self.ligands[ligand].load_crystal_pose()
            else:
                self.ligands[ligand].load_poses(load_fp)
            self.num_poses[ligand] = len(self.ligands[ligand].poses)

class LigandManager:
    """
    Manages docking results, fps, and MCSSs for all ligands associated
    with a given protein.

    This class defines path to protein directory.

    Parameterized by (protein).
    """
    def __init__(self, protein, root, struct = 'First'):
        """
        root (string): Path to protein's data directory.
        """
        self.protein = protein
        self.root = root
        
        # Initialize ligand info.
        self.chembl_info = load_chembl_proc(self.root)
        self.u_ligs, self.dup_ligs = read_duplicates(self.root)
        self.all_ligs = self.prepped()
        self.pdb = self.unique(sorted([l for l in self.all_ligs
                                      if l[:6] != 'CHEMBL'],
                                      reverse = struct == 'Last'))

        # Set default structure.
        self.st = None
        if not os.path.exists('{}/docking/grids'.format(self.root)): return
        self.grids = sorted([l for l in os.listdir('{}/docking/grids'.format(self.root)) if l[0] != '.'])
        if not self.grids: return
        
        if struct is 'First':
            self.st = self.grids[0] 
        elif struct is 'Last':
            self.st = self.grids[-1]
        else:
            assert False

        if self.protein == 'AR':
            self.st = '2AXA'

        if self.protein == 'NR3C1':
            self.st = '3BQD'

        if self.protein == 'NR3C2':
            self.st = '3WFF'

        self.mcss = MCSSController(self)
        self.helpers = {}

    def get_xdocked_ligands(self, num):
        ligands = self.docked(self.pdb)[:num+1]
        self_docked = self.st+'_lig'
        if self_docked in ligands:
           ligands.remove(self_docked)
        elif ligands:
            ligands.pop(-1)
        return ligands

    def docked(self, ligands, st=None):
        if st == None: st = self.st
        return [ligand for ligand in ligands
                if os.path.exists("{}/{}-to-{}_pv.maegz".format(
                                  Ligand(ligand, Docking(self.root, st).dock_dir,
                                         '', st).glide_path,
                                  ligand, st
                                  ))]

    def prepped(self):
        ligdir = '{}/ligands/prepared_ligands'.format(self.root)
        if not os.path.exists(ligdir): return set([])
        return set([l for l in os.listdir(ligdir) if os.path.exists('{}/{}/{}_out.mae'.format(ligdir,l,l))])

    def chembl(self):
        def valid(l):
            filters = [
                lambda x,ci: ci[x].ki is not None and ci[x].ki <= 1000,
                lambda x,ci: ci[x].mw is not None and ci[x].mw <= 800,
                lambda x,ci: ci[x].macrocycle is not None and not ci[x].macrocycle,
                lambda x,ci: ci[x].valid_stereo is not None and ci[x].valid_stereo
            ]
            return all(f(l, self.chembl_info) for f in filters)
        return [l for l in self.all_ligs if l in self.chembl_info and valid(l)]

    def unique(self, l_list):
        """
        Removes duplicates from l_list.
        
        If identical ligands are found, the one that appears first in l_list will be kept.
        """
        if not self.u_ligs:
            print('duplicates not loaded')
            return l_list
        unique_ligs = []
        exclude = set([])
        for l in l_list:
            if l in exclude: continue
            unique_ligs.append(l)
            if l in self.u_ligs: continue
            for l2, dups in self.dup_ligs.items():
                if l in dups:
                    exclude.update(dups)
                    break
            else:
                print('uh oh, ligand not found in unique or duplicates...', l)
        return unique_ligs

    def pick_helpers(self, maxnum=20, num_chembl=20):
        parent = 'chembl/helpers'
        os.system('mkdir -p {}'.format(parent))
        ki_sort = lambda c: self.chembl_info[c].ki

        # filters
        all_options = [
            'best_affinity.txt',
            'best_mcss.txt',
            'best_affinity_diverse.txt'
        ]

        for f in all_options:
            fpath = '{}/{}'.format(parent, f)
            if not os.path.exists(fpath):
                print('picking chembl ligands', f)
                chembl_ligs = sorted(self.chembl())
                with open(fpath, 'w') as fi:
                    for q in self.get_xdocked_ligands(maxnum):
                        print(q)
                        if 'mcss' in f:
                            self.mcss.load_mcss()
                            sorted_helpers = sorted(chembl_ligs, key=ki_sort)
                            sorted_helpers = self.mcss.sort_by_mcss(q, sorted_helpers)
                        elif 'affinity' in f:
                            sorted_helpers = sorted(chembl_ligs, key=ki_sort)

                        unique = self.unique([q]+sorted_helpers)

                        if f == 'best_affinity_diverse.txt':
                            self.mcss.load_mcss()
                            _unique = []
                            for helper in unique:
                                for _helper in _unique:
                                    not_query = _helper != q
                                    sim = self.mcss.get_mcss_size(helper, _helper,compute=True) > 0.8
                                    if sim and not_query:
                                        break
                                else:
                                    _unique += [helper]

                                print(len(_unique))
                                if len(_unique) == num_chembl+1:
                                    break
                            unique = _unique
                        fi.write('{}:{}\n'.format(q, ','.join(unique[1:num_chembl+1])))

    def load_helpers(self):
        fpath = 'chembl/helpers'
        if dirpath is not None:
            fpath = '{}/{}'.format(dirpath, fpath)

        helpers = {}
        for fname in os.listdir(fpath):
            if fname[0] == '.' or fname.split('.')[-1] != 'txt': continue
            helpers[fname] = {}
            with open('{}/{}'.format(fpath, fname)) as f:
                for line in f:
                    q, chembl = line.strip().split(':')
                    helpers[fname][q] = chembl.split(',')
        return helpers

    def get_grids(self):
        return self.grids

class Protein:
    """
    Collection of ligands and docking results for a given protein.
    """
    def __init__(self, protein, struct = 'First'):
        self.root = "{}/{}".format(shared_paths['data'], protein)
        self.lm = LigandManager(protein, self.root, struct)
        # Useful to be able to reference this before loading data.
        self.docking = {self.lm.st: Docking(self.root, self.lm.st)}

    def load_docking(self, ligands, load_fp=False, load_crystal = False,
                     load_mcss = False, st = None):
        if load_crystal:
            assert all('crystal' in ligand for ligand in ligands)
        else:
            assert not any('crystal' in ligand for ligand in ligands)

        if st == None: st = self.lm.st

        if st not in self.docking:
            self.docking[st] = Docking(self.root, st)
        self.docking[st].load(ligands, load_fp, load_crystal)

        if load_mcss:
            self.lm.mcss.load_rmsds(ligands+list(self.docking[st].ligands.keys()),
                                    shared_paths['stats']['max_poses'])
