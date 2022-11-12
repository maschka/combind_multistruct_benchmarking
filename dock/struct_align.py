import os
from subprocess import run

command = ('$SCHRODINGER/utilities/structalign '
           '-asl "(not chain. L and not atom.element H) and (fillres within {0} chain. L)" '
           '-asl_mobile "(not chain. L and not atom.element H) and (fillres within {0} chain. L)" '
           '{1} {2}')

#hack to not include ligands here!
def should_include_lig(struct):
    if "_" in struct:
        print('modeled_struct!')
        return True
    return True

def align_successful(out_dir, struct):

    if not os.path.exists('{}/{}/rot-{}_query.mae'.format(out_dir, struct, struct)):
        return False
    
    if os.path.exists('{}/{}/{}_template.mae'.format(out_dir, struct, struct)):
        return True # query = template so we don't need to check alignment

    with open('{}/{}/align.out'.format(out_dir, struct), 'r') as f:
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

def align_successful_all(out_dir, struct, template):

    if not os.path.exists('{}/{}/rot-{}_query_to_{}.mae'.format(out_dir, struct, struct,template)):
        return False
    
    # if os.path.exists('{}/{}/{}_template.mae'.format(out_dir, struct, struct)):
        # return True # query = template so we don't need to check alignment

    with open('{}/{}/align_to_{}.out'.format(out_dir, struct,template), 'r') as f:
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

def struct_align(template, structs, dist=15.0, retry=True,
                 processed_out='structures/processed/{pdb}/{pdb}_out.mae',
                 align_dir='structures/aligned'):
    
    template_path = processed_out.format(pdb=template)
    if not os.path.exists(template_path):
        print('template not processed', template_path)
        return

    for struct in structs:
        query_path = processed_out.format(pdb=struct)
        if not os.path.exists(query_path) or align_successful(align_dir, struct):
            continue

        print('align', struct, template)

        os.system('mkdir -p {}'.format(align_dir))
        os.system('rm -rf {}/{}'.format(align_dir, struct))
        os.system('mkdir -p {}/{}'.format(align_dir, struct))

        _workdir = '{}/{}'.format(align_dir, struct)
        _template_fname = '{}_template.mae'.format(template)
        _query_fname = '{}_query.mae'.format(struct, template)

        os.system('cp {} {}/{}'.format(template_path, _workdir, _template_fname))
        os.system('cp {} {}/{}'.format(query_path, _workdir, _query_fname))

        with open('{}/align_in.sh'.format(_workdir), 'w') as f:
            f.write(command.format(dist, _template_fname, _query_fname))
        run('sh align_in.sh > align.out', shell=True, cwd=_workdir)

        if retry and not align_successful(align_dir, struct):
            print('Alignment failed. Trying again with a larger radius.')
            struct_align(template, [struct], dist=25.0, retry=False,
                         processed_out=processed_out, align_dir=align_dir)



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
        if align_successful_all(align_dir, struct,template):
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

        if retry and not align_successful_all(align_dir, struct,template):
            print('Alignment failed. Trying again with a larger radius.')
            struct_align_all(template, [struct], dist=25.0, retry=False,
                         processed_out=processed_out, align_dir=align_dir)
