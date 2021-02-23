from multiprocessing import Pool
import os
import numpy as np
from schrodinger.structure import StructureReader

def np_load(fname, halt=True, delete=False):
    fname = os.path.abspath(fname)
    try:
        return np.load(fname)
    except ValueError as e:
        m = 'Cannot load file containing pickled data when allow_pickle=False'
        if m in str(e):
            print('{} is corrupt. Regenerate and try again.'.format(fname))
            if delete:
                os.remove(fname)
        else:
            print("Can't open {}".format(fname))
            print(str(e))

        if halt:
            exit()

def pv_path(root, name):
    if '_native' in name:
        name = name.replace('_native', '')
        return '{}/{}/{}_native_pv.maegz'.format(root, name, name)
    return '{}/{}/{}_pv.maegz'.format(root, name, name)

def get_pose(pv, pose):
    with StructureReader(pv) as sts:
        for _ in range(pose+1):
            next(sts)
        st = next(sts)
    return st

def basename(path):
    x = os.path.basename(path)
    x = os.path.splitext(x)[0]
    return x

def mp(function, unfinished, processes):
    if unfinished:
        with Pool(processes=processes) as pool:
            x = pool.starmap(function, unfinished)
        return x

def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)
