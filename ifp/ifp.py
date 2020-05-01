import click
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MaeMolSupplier
import gzip
import rdkit

################################################################################

def read_mae(mae, poses):
    if mae[-2:] == 'gz':
        with gzip.open(mae) as fp:
            mols = []
            for mol in MaeMolSupplier(fp, removeHs=False):
                if mol is None:
                    print(mae)
                    exit()
                mols += [mol]
    else:
        mols = [mol for mol in MaeMolSupplier(mae, removeHs=False)]
    return mols

def resname(atom):
    info = atom.GetPDBResidueInfo()
    return ':'.join(map(lambda x: str(x).strip(),
                        [info.GetChainId(), str(info.GetResidueNumber()),
                         info.GetResidueName(), info.GetInsertionCode()]))

def atomname(atom):
    pdb = atom.GetPDBResidueInfo()
    if pdb is None:
        return str(atom.GetIdx())
    return pdb.GetName().strip()

def coords(atom):
    return atom.GetOwningMol().GetConformer(0).GetAtomPosition(atom.GetIdx())

def distance(atom1, atom2):
    return coords(atom1).Distance(coords(atom2))

def angle(atom1, atom2, atom3):
    v1 = coords(atom1) - coords(atom2)
    v3 = coords(atom3) - coords(atom2)
    return v1.AngleTo(v3) * 180.0 / np.pi

################################################################################

def _get_bonded_hydrogens(atom):
    hydrogens = []
    for bond in atom.GetBonds():
        if bond.GetBeginAtomIdx() != atom.GetIdx():
            hydrogen = bond.GetBeginAtom()
        else:
            hydrogen = bond.GetEndAtom()
            
        if hydrogen.GetAtomicNum() == 1:
            hydrogens += [hydrogen]
    return hydrogens

def _is_donor(atom):
    if atom.GetAtomicNum() in [7, 8]:
        if _get_bonded_hydrogens(atom):
            return True
    return False

def _is_acceptor(atom):
    if atom.GetAtomicNum() == 8:
        return True
    if atom.GetAtomicNum() == 7 and atom.GetExplicitValence() == 3:
        return True
    return False

def _hbond_hydrogen_angle(acceptor, donor):
    best_angle, best_hydrogen = 0, None
    for hydrogen in _get_bonded_hydrogens(donor):
        _angle = angle(donor, hydrogen, acceptor)
        if _angle > best_angle:
            best_angle = _angle
            best_hydrogen = hydrogen
    return best_hydrogen, best_angle

def _hbond_compute(donor_mol, acceptor_mol, settings, protein_is_donor):
    donors = [atom for atom in donor_mol if _is_donor(atom)]
    acceptors = [atom for atom in acceptor_mol if _is_acceptor(atom)]

    hbonds = []
    for donor in donors:
        for acceptor in acceptors:
            dist = distance(acceptor, donor)
            hydrogen, angle = _hbond_hydrogen_angle(acceptor, donor)
            if (dist < settings['hbond_dist_cut']
                and angle > settings['hbond_angle_cut']):

                if protein_is_donor:
                    label = 'hbond_donor'
                    protein_atom = donor
                    ligand_atom = acceptor
                else:
                    label = 'hbond_acceptor'
                    protein_atom = acceptor
                    ligand_atom = donor

                hbonds += [{'label': label,
                            'protein_res': resname(protein_atom),
                            'protein_atom': atomname(protein_atom),
                            'ligand_atom': atomname(ligand_atom),
                            'dist': dist,
                            'angle': angle,
                            'hydrogen': atomname(hydrogen)}]
    return hbonds

def hbond_compute(protein, ligand, settings):
    donor = _hbond_compute(protein, ligand, settings, True)
    acceptor = _hbond_compute(ligand, protein, settings, False)
    return acceptor + donor

################################################################################

def _symetric_charged_protein_atoms(protein):
    protein_groups = {}
    for protein_atom in protein:
        if atomname(protein_atom) in ['OD1', 'OD2', 'OE1', 'OE2', 'NH1', 'NH2']:
            if resname(protein_atom) not in protein_groups:
                protein_groups[resname(protein_atom)] = []
            protein_groups[(resname(protein_atom))] += [protein_atom]
    return protein_groups

def _symetric_charged_ligand_atoms(ligand):
    ligand_groups = {}
    smartss = [('[CX3](=O)[O-]', 2, [1, 2]),
               ('[CX3](=[NH2X3+])[NH2X3]', 1, [1, 2])]

    idx_to_atom = {atom.GetIdx(): atom for atom in ligand}

    for smarts, k, v in smartss:
        mol = Chem.MolFromSmarts(smarts)
        matches = ligand[0].GetOwningMol().GetSubstructMatches(mol)
        for match in matches:
            ligand_groups[match[k]] = [idx_to_atom[match[_v]] for _v in v]
    return ligand_groups

def saltbridge_compute(protein, ligand, settings):
    # Note that much of the complexity here stems from taking into account
    # symetric atoms. Specifically for carboxylate and guanidinium groups,
    # we consider not just the atom that is arbitrarily assigned a formal
    # charge, but also the atom that is charged in the other resonance
    # structure.
    protein_groups = _symetric_charged_protein_atoms(protein)
    ligand_groups = _symetric_charged_ligand_atoms(ligand)

    saltbridges = []
    for protein_atom in protein:
        if protein_atom.GetFormalCharge() == 0: continue
        for ligand_atom in ligand:
            if ligand_atom.GetFormalCharge() == 0: continue
            lig_charge = ligand_atom.GetFormalCharge()
            protein_charge = protein_atom.GetFormalCharge()
            if lig_charge * protein_charge >= 0: continue

            # Expand protein_atom and ligand_atom to all symetric atoms
            # ... think carboxylates and guanidiniums.
            if ligand_atom.GetIdx() in ligand_groups:
                ligand_atoms = ligand_groups[ligand_atom.GetIdx()]
            else:
                ligand_atoms = [ligand_atom]
            if resname(protein_atom) in protein_groups:
                protein_atoms = protein_groups[resname(protein_atom)]
            else:
                protein_atoms = [protein_atom]
            
            # Get minimum distance between any pair of protein and ligand
            # atoms in the groups.
            dist = float('inf')
            for _ligand_atom in ligand_atoms:
                for _protein_atom in protein_atoms:
                    _dist = distance(_protein_atom, _ligand_atom)
                    if _dist < dist:
                        dist = _dist
                        protein_atom = _protein_atom
                        ligand_atom = _ligand_atom

            if dist < settings['sb_dist_cut']:
                saltbridges += [{'label': 'saltbridge',
                                 'protein_res': resname(protein_atom),
                                 'protein_atom': atomname(protein_atom),
                                 'ligand_atom': atomname(ligand_atom),
                                 'dist': dist}]
    return saltbridges

################################################################################

def contact_compute(protein, ligand, settings):
    contacts = []
    for protein_atom in protein:
        if protein_atom.GetAtomicNum() not in settings['nonpolar']: continue
        for ligand_atom in ligand:
            if ligand_atom.GetAtomicNum() not in settings['nonpolar']: continue
    
            vdw = settings['nonpolar'][ligand_atom.GetAtomicNum()]
            vdw += settings['nonpolar'][protein_atom.GetAtomicNum()]
            dist = distance(protein_atom, ligand_atom)

            if dist < vdw*settings['contact_scale_cut']:
                contacts += [{'label': 'contact',
                              'protein_res': resname(protein_atom),
                              'protein_atom': atomname(protein_atom),
                              'ligand_atom': atomname(ligand_atom),
                              'dist': dist,
                              'vdw': vdw}]
    return contacts

################################################################################

def _protein_coords(protein):
    return np.array([coords(atom) for atom in protein])

def _get_relevent_atoms(protein, ligand, protein_coords, overall_cut):
    ligand_coords  = np.array([coords(atom) for atom in ligand])
    
    protein_coords = np.expand_dims(protein_coords, 1)
    ligand_coords = np.expand_dims(ligand_coords, 0)
    
    displacements = ligand_coords - protein_coords
    distances = np.linalg.norm(displacements, axis=2)
    distances = distances.min(axis=1)
    return protein[distances < overall_cut]

def _fingerprint(protein, ligand, protein_coords, settings):
    # Here protein and ligand should be a numpy array of rdkit atoms
    # and protein_coords should be a numpy array of protein atomic
    # coordinates.
    protein = _get_relevent_atoms(protein, ligand, protein_coords, settings['overall_cut'])

    fp  = hbond_compute(protein, ligand, settings)
    fp += saltbridge_compute(protein, ligand, settings)
    fp += contact_compute(protein, ligand, settings)
    return pd.DataFrame.from_dict(fp)

def fingerprint(protein, ligand, settings):
    _ligand = np.array([atom for atom in ligand.GetAtoms() if atom.GetAtomicNum() != 1])
    _protein = np.array([atom for atom in protein.GetAtoms() if atom.GetAtomicNum() != 1])
    protein_coords = _protein_coords(_protein)
    return _fingerprint(_protein, _ligand, protein_coords, settings)

def fingerprint_poseviewer(protein, ligands, settings):

    # Get heavy atoms and position matrix to optimize look-up of relevent atoms.
    _ligands = [np.array([atom for atom in ligand.GetAtoms() if atom.GetAtomicNum() != 1])
                for ligand in ligands]
    _protein = np.array([atom for atom in protein.GetAtoms() if atom.GetAtomicNum() != 1])
    protein_coords = _protein_coords(_protein)

    fps = []
    for i, ligand in enumerate(_ligands):
        fps += [_fingerprint(_protein, ligand, protein_coords, settings)]
        fps[-1]['pose'] = i
    fps = pd.concat(fps, ignore_index=True, sort=False)
    if 'hydrogen' not in fps:
        fps['hydrogen'] = ''
    fps.loc[fps['hydrogen'].isna(), 'hydrogen'] = ''
    return fps

################################################################################

def _piecewise(data, opt, cut):
    slope = 1 / (cut-opt)
    intercept = cut * slope

    data = intercept - slope * data
    data[data > 1] = 1
    data[data < 0] = 0
    return data

def _groupby_subset(df, index, col):
    return df[index+[col]].groupby(index)

def compute_scores(raw, settings):
    scores = []
    for label, group in raw.groupby('label'):
        group = group.copy()
        if label == 'contact':
            group['score'] = _piecewise(group['dist'] / group['vdw'],
                                        settings['contact_scale_opt'],
                                        settings['contact_scale_cut'])
        elif label == 'saltbridge':
            group['score'] = _piecewise(group['dist'],
                                        settings['sb_dist_opt'],
                                        settings['sb_dist_cut'])
        elif label in ['hbond_donor', 'hbond_acceptor']:
            group['score'] = (  _piecewise(group['dist'],
                                           settings['hbond_dist_opt'],
                                           settings['hbond_dist_cut'])
                              * _piecewise(180 - group['angle'],
                                           settings['hbond_angle_opt'],
                                           settings['hbond_angle_cut']))

            # One hydrogen bond per hydrogen
            if label == 'hbond_donor':
                idx = _groupby_subset(group,
                                      ['pose', 'protein_res', 'hydrogen'],
                                      'score').idxmax()
            else:
                idx = _groupby_subset(group,
                                      ['pose', 'hydrogen'],
                                      'score').idxmax()
            idx = idx['score']
            group = group.loc[idx]
        group = _groupby_subset(group,
                                ['pose', 'label', 'protein_res'],
                                'score').sum()
        scores += [group]
    return pd.concat(scores).sort_index()

################################################################################

def ifp(settings, input_file, output_file, poses):
    settings['nonpolar'] = {6:1.7, 9:1.47, 17:1.75, 35:1.85, 53:1.98}
    settings['overall_cut'] = max(settings['hbond_dist_cut'], settings['sb_dist_cut'],
                                  settings['contact_scale_cut']*2*max(settings['nonpolar'].values()))

    protein, *ligands = read_mae(input_file, poses)
    ligands = ligands[:poses]

    fps = fingerprint_poseviewer(protein, ligands, settings)
    scores = compute_scores(fps, settings)

    fps = fps.set_index(['pose', 'label', 'protein_res', 'protein_atom', 'ligand_atom'])
    fps = fps.sort_index()
    base = output_file.split('.')
    base, ext = base[:-1], base[-1]
    raw_file = '.'.join(base) + '_raw.' + ext

    fps.to_csv(raw_file)
    scores.to_csv(output_file)
