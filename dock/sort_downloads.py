import os
import sys
import pandas as pd
from glob import glob

from datetime import date

from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.structutils.analyze import evaluate_asl

# taken from schrodinger "solvent" asl definition
solvent = ['hoh', 'spc', 't4p', 't3p', 't5p', 't4pe', 'dod', 'GOL', 'EDO', 'DMS', 'PEG', 'MPD', 
           'PG4', 'BME', 'PGE', 'IPA', '1PE', 'EOH', 'P6G', 'DMF', 'MOH', 'POL', '2PE', 'PE4',
           'ETA', '15P', '12P', 'P33', 'PE5', 'TBU', 'CCN', 'CXE', '7PE', 'NME', 'PE3', 'PE6', 'NHE',
           'ACN', 'P4C', 'PE8', 'NEH', 'XPE', '1BO', 'N8E', 'DMN', 'CE1', 'PYE', 'C10', 'SBT', 'MB3']
solvent = [x.lower() for x in solvent]

ignore_ligands = solvent + ['', 'dtt', 'epe', 'mes', 'tla', 'tar', 'nmm', 
                 'eu','au','edt', 'ola','olb','olc','css', 'cso', 'clr', 'pg0', 'srt', 'cxs', 'tpo', 'sog',
                 'ccs','cme','ben','glc','mli','ocs','aly','sgm','ptr','csd','kcx','ste','cit','iod','hto',
                 'hez','jzr','cps','bog','acp','anp','adp','atp','ags', 'glu']

ignore_titles = ['mutant', 'mutation', 't877a', 'w741l', 'cryptic']#, 'fragment']#, 'allosteric']

ignore_methods = ['solution nmr']


# These should be specified on the fly in pdb_procesing_manual.txt
# orthosteric = ['dht','tes'] # ar       
# allosteric = ['8vb','8z5','z24','imw', # plk1
#               '2an', # cdk2
#               'ld2'] # mr
# ignore_pdb = ['1h00','1h01','1h07','1h08', # cdk2 2 small molecules on top of each other?
#               '5tlx','5kcf','5kct','5kra','5tlp', # era multiple small molecules per structure
#               '4h71'] # plk1 allosteric site
    
def get_ligs(structs, log):
    log.write('Begin dismabiguating ligands.\n')
    ligs = {}
    for pdbid, pdblist in sorted(structs.items(),key=lambda x:x[0]):
        if len(set([p.lig for p in pdblist])) != 1:
            log.write("Ignoring PDB {}: Unable to choose ligand from {}.\n".format(pdbid, set([p.lig for p in pdblist])))
            continue
        if pdblist[0].lig not in ligs: ligs[pdblist[0].lig] = []    
        first_chain = sorted(pdblist, key=lambda x:x.chain)[0] 
        ligs[pdblist[0].lig].append(first_chain)
    log.write('End disambiguating ligands')
    return ligs

def write_files(ligs, log):
    log.write('Begin file writing.\n')
    for lname, l_copies in sorted(ligs.items(),key=lambda x:x[1][0].pdb):
        if any(os.path.exists('structures/raw_files/{}_lig.mae'.format(l.pdb.upper())) for l in l_copies):
            log.write("Ignoring PDBs {}: Ligand already present in raw_files\n".format(','.join(l_copies)))

        l_copies.sort(key=lambda x: x.date)

        # find the first non-covalent copy
        for l in l_copies: 
            prot_st = next(StructureReader('structures/downloads/{}.pdb'.format(l.pdb)))
            lig_st = None
            for res in prot_st.residue:
                if res.pdbres.strip().lower() == lname and res.chain.strip().lower() == l.chain:
                    molnum = res.molecule_number
                    if len([m.residue for m in prot_st.molecule if m.number == molnum][0]) == 1:
                        lig_st = prot_st.extract([a.index for a in res.atom])
                    break
            if lig_st is None:
                log.write("Ignoring PDB {}: Ligand covalently bound\n".format(l.pdb))
                continue
            
            for pdb in l_copies[l_copies.index(l)+1:]:
                log.write("Ignoring PDB {}: Ligand redundant with {}\n".format(pdb.pdb, l.pdb))

            log.write("Writing PDB {}.\n".format(l.pdb))
            to_delete = evaluate_asl(prot_st, 'solvent or res.pt {}'.format(lname.upper()))
            prot_st.deleteAtoms(to_delete)
                
            prot_wr = StructureWriter('structures/raw_files/{}_prot.mae'.format(l.pdb.upper()))
            prot_wr.append(prot_st)
            prot_wr.close()
            
            lig_wr = StructureWriter('structures/raw_files/{}_lig.mae'.format(l.pdb.upper()))
            lig_wr.append(lig_st)
            lig_wr.close()
            break
    log.write('End file writing.\n')


class Chain:
    """
    Specifies ligands associated with each chain of a PDB structure.
    """
    template = "{}:{}"
    def __init__(self, chain_id, ligands):
        self.chain_id = chain_id
        self.ligands = ligands

    def add_ligands(self, ligands):
        self.ligands += ligands

    def dumps(self):
        return self.template.format(self.chain_id, ','.join(self.ligands))

    @classmethod
    def loads(cls, S):
        chain_id, ligands = S.split(':')
        return Chain(chain_id, ligands.split(','))

class PDB:
    """
    Specifications whether to use a protein and, if so, which
    ligand and chain to use.
    """
    template = "pdb={}\tdate={}\tchains={}\tdelete={}\ttitle={}\tcomments={}"

    def __init__(self, pdb, title, date):
        """
        pdb: 4 letter PDB code.
        chains: dict mapping chains to associated ligands
        """
        self.pdb = pdb
        self.title = title
        self.date = date
        
        self.chains = {}
        self.delete = []
        self.comments = {}

    def add_chain(self, chain_id, ligand):
        if chain_id not in self.chains: self.chains[chain_id] = Chain(chain_id, [])
        self.chains[chain_id].add_ligands([ligand])

    def add_comment(self,  chain_id, message):
        if chain_id not in self.comments: self.comments[chain_id] = []
        self.comments[chain_id] += [message]

    def ligands(self):
        return [(chain.chain_id, ligand) for chain in self.chains.values() for ligand in chain.ligands if ligand]

    def dumps(self):
        return self.template.format(self.pdb, self.date,
                               ';'.join(chain.dumps() for chain in self.chains.values()),
                               ','.join(self.delete), self.title,
                               ';'.join("{}:{}".format(chain, ','.join(comments))
                                        for chain, comments in self.comments.items())
                                )
    @classmethod
    def loads(cls, S):
        for s in S.strip().split('\t'):
            k, v = s.split('=')
            if k == 'pdb': pdb = v
            if k == 'date': date = v
            if k == 'chains': chains = [Chain.loads(chain) for chain in v.split(';')] if v else []
            if k == 'delete': delete = v.split(',')
            if k == 'title': title = v
            if k == 'comments': comments = {chain.split(':')[0]: chain.split(':')[1].split(',') for chain in v.split(';')} if v else {}
        pdb = PDB(pdb, title, date)
        pdb.chains = {chain.chain_id:chain for chain in chains}
        pdb.delete = delete
        pdb.comments = comments
        return pdb

def read_initial_pdbs():
    # load in csv file
    pdbs = {}
    for csv in glob('structures/downloads/*.csv'):
        uniprot = csv.split('/')[-1].split('.')[0].lower()

        df = pd.read_csv(csv)
        for col in df:
            for col in df:
                if type(df[col][0]) == str:
                    df[col] = df[col].str.lower()

        for i, line in df.iterrows():
            m,d,y = line['Rel. Date'].split('/')
            st_date = date(int(y),int(m),int(d))

            if line['PDB ID'] not in pdbs:
                pdbs[line['PDB ID']] = PDB(line['PDB ID'], line['Structure Title'], st_date)

            if line['Exp. Method'] in ignore_methods:
                pdbs[line['PDB ID']].add_comment(line['Chain ID'], "Method {} in ignore_methods".format(line['Exp. Method']))
                continue
            if uniprot not in line['Uniprot Acc'].split(','):
                pdbs[line['PDB ID']].add_comment(line['Chain ID'], "Uniprot not target uniprot")
                continue
            if any(x in line['Structure Title'] for x in ignore_titles):
                pdbs[line['PDB ID']].add_comment(line['Chain ID'], "Title contains one of ignore_titles".format(line['Structure Title']))
                continue
            if line['Ligand ID'] in ignore_ligands:
                pdbs[line['PDB ID']].add_comment(line['Chain ID'], "Ligand {} in ignore_ligands".format(line['Ligand ID']))
                continue

            if line['Ligand ID'] in [ligand for chain, ligand in pdbs[line['PDB ID']].ligands()]:
                pdbs[line['PDB ID']].add_comment(line[chain_i], "Ligand {} present in earlier chain".format(line['Ligand ID']))
                continue

            if not (100 < float(line['Ligand MW']) < 1000):
                pdbs[line['PDB ID']].add_comment(line['Chain ID'], "Ligand {} of MW {} not right size".format(line['Ligand ID'], line['Ligand MW']))
                continue
                
            pdbs[line['PDB ID']].add_chain(line['Chain ID'], line['Ligand ID'])
                
    return pdbs

def mark_covalently_bound(pdbs):
    for pdb_id, pdb in pdbs.items():
        prot_st = next(StructureReader('structures/downloads/{}.pdb'.format(pdb.pdb)))
        for chain_id, ligand in pdb.ligands():
            found = False
            for res in prot_st.residue:
                if res.pdbres.strip().lower() == ligand and res.chain.strip().lower() == chain_id:
                    molnum = res.molecule_number
                    if len([m.residue for m in prot_st.molecule if m.number == molnum][0]) == 1:
                        found = True
            if not found:
                pdb.chains[chain_id].ligands.remove(ligand)
                pdb.add_comment(chain_id, "Ligand {} covalently bound".format(ligand))

def mark_duplicates(pdbs):
    """
    Currently ignoring the case where
    """
    ligands = {}
    for pdb_id, pdb in pdbs.items():
        for chain_id, ligand in pdb.ligands():
            if ligand not in ligands: ligands[ligand] = []
            ligands[ligand] += [pdb_id]

    for ligand in ligands:
        ligands[ligand] = sorted(ligands[ligand], key=lambda x: pdbs[x].date)

    # Check if each structure has any unique ligands
    for pdb_id, pdb in pdbs.items():
        unique = False
        for chain_id, ligand in pdb.ligands():
            if len(ligands[ligand]) == 1 or ligands[ligand][0] == pdb_id:
                unique = True
            else:
                pdb.add_comment(chain_id, "Ligand {} is redundant with PDB {}".format(ligand, ligands[ligand][0]))
        # Remove all ligands if none are unique
        if not unique:
            for chain in pdb.chains.values():
                chain.ligands = []

def write_pdbs(pdbs, path):
    """
    Write structures to PATH for each pdb with 1+ ligands.
    If a structure has multiple ligands, write multiple entries
    with basename PDB_LIG instead of just PDB.
    """
    for pdb_name, pdb in pdbs.items():
        for chain_id, ligand in pdb.ligands():
            prot_st = next(StructureReader('structures/downloads/{}.pdb'.format(pdb.pdb)))
            lig_st = None
            for res in prot_st.residue:
                if res.pdbres.strip().lower() == ligand and res.chain.strip().lower() == chain_id:
                    molnum = res.molecule_number
                    if len([m.residue for m in prot_st.molecule if m.number == molnum][0]) == 1:
                        lig_st = prot_st.extract([a.index for a in res.atom])

            # TODO: Add pdb.to_delete to this selection
            to_delete = evaluate_asl(prot_st, 'solvent or res.pt {}'.format(ligand.upper()))
            prot_st.deleteAtoms(to_delete)

            basename = pdb.pdb.upper() if len(pdb.ligands()) == 1 else "{}_{}".format(pdb.pdb.upper(), ligand.upper())
                
            prot_wr = StructureWriter('structures/{}/{}_prot.mae'.format(path, basename))
            prot_wr.append(prot_st)
            prot_wr.close()
            
            lig_wr = StructureWriter('structures/{}/{}_lig.mae'.format(path, basename))
            lig_wr.append(lig_st)
            lig_wr.close()

def write_log(pdbs):
    with open('structures/initial_pdb_processing.txt', 'w') as fp:
        # ambiguous
        for pdb_id, pdb in pdbs.items():
            if len(pdb.ligands()) > 1:
                fp.write(pdb.dumps()+'\n')
        # unambiguous
        for pdb_id, pdb in pdbs.items():
            if len(pdb.ligands()) == 1:
                fp.write(pdb.dumps()+'\n')
        # no ligands
        for pdb_id, pdb in pdbs.items():
            if len(pdb.ligands()) < 1:
                fp.write(pdb.dumps()+'\n')

def read_log():
    pdbs = {}
    with open('structures/final_pdb_processing.txt') as fp:
        for line in fp:
            pdb = PDB.loads(line)
            pdbs[pdb] = pdb
    return pdbs

def check_unambiguous(pdbs):
    for pdb_id, pdb in pdbs.items():
        assert len(pdb.ligands()) <= 1, pdb.pdb

    ligands = set()
    for pdb_id, pdb in pdbs.items():
        for chain, ligand in pdb.ligands():
            assert ligand not in ligands, ligand
            ligands.add(ligand)

def initial_processing():
    """
    Output file initial_pdb_processing.txt containing a line
    for each pdb file.
    """
    pdbs = read_initial_pdbs()
    mark_covalently_bound(pdbs)
    mark_duplicates(pdbs)
    write_pdbs(pdbs, 'initial_files')
    write_log(pdbs)

def final_processing():
    """
    Reads the file final_pdb_processing.txt, which is expected to
    have unambiguous instructions on what pdb structures, ligands to use,
    and het atoms to delete outside of the ones specified above.
    """
    pdbs = read_log()
    check_unambiguous(pdbs)
    write_pdbs(pdbs, 'raw_files')


def sort_downloads():
    """
    Workflow is as follows:
    (1) Download pdb files + report describing thier contents
    (2) Run initial processing, which attempts to determine the
        orthosteric ligand based on simple rules about known names
        of solvent molecules and the ligand molecular weights.
        Writes all possible options to initial_files/*
    (3) Resolve remaining ambiguities manually by inspecting the
        file initial_pdb_processing.txt and deleting entries for
        ligands that should not be considered. Write this updated
        file to final_pdb_processing.pdb.
    (4) Run final processing, which writes the final structures
        and ligands to raw_files/*
    """
    if os.path.exists('structures/final_pdb_processing.txt'):
        print('Running final processing')
        final_processing()
    else:
        print('Running initial processing')
        initial_processing()
