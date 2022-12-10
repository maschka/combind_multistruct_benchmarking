import os
from schrodinger.structure import StructureReader
import shutil

#hack to not include ligands here!
def should_include_lig(struct):
    if "_" in struct:
        print('modeled_struct!')
        return True
    return True

def split_complex(st, pdb_id):
    os.system('mkdir -p structures/proteins structures/ligands')
    lig_path = 'structures/ligands/{}_lig.mae'.format(pdb_id)
    prot_path = 'structures/proteins/{}_prot.mae'.format(pdb_id)
    print('split?')
    #need to copy from raw, because wasn't minimized! I'm also copyin unrotated structure here!
   
    # if not os.path.exists(lig_path) and len([a.index for a in st.atom if a.chain == 'L']) > 0:
    if not should_include_lig(pdb_id):
        # print('here?')
        if os.path.exists('structures/raw/{}_lig.mae'.format(pdb_id)):
            shutil.copy( 'structures/raw/{}_lig.mae'.format(pdb_id), 'structures/ligands/{}_lig.mae'.format(pdb_id))
    else: 
        lig_st = st.extract([a.index for a in st.atom if a.chain == 'L'])
        lig_st.title = '{}_lig'.format(pdb_id)
        lig_st.write(lig_path)
    
    # if not os.path.exists(prot_path):
    prot_st = st.extract([a.index for a in st.atom if a.chain != 'L'])
    prot_st.title = '{}_prot'.format(pdb_id)
    prot_st.write(prot_path)

#version used rotation, but I don't want that!
# def struct_sort(structs):
#     for struct in structs:
#         opt_complex = 'structures/aligned/{}/rot-{}_query.mae'.format(struct, struct)

#         if os.path.exists(opt_complex):
#             comp_st = next(StructureReader(opt_complex))
#             split_complex(comp_st, struct)

def struct_sort(structs):
    for struct in structs:
        #using unrotated query here. makes it easier to align everythin later.
	# john 11/13 changed to merged.mae for test 
        # opt_complex = 'structures/aligned/{}/{}_query.mae'.format(struct, struct)
        opt_complex = 'structures/processed/{}/{}_out.mae'.format(struct, struct)
        if os.path.exists(opt_complex):
            print('help')
            comp_st = next(StructureReader(opt_complex))
            split_complex(comp_st, struct)
