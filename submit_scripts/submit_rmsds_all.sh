#!/bin/bash
#for struct in ../data/*/structures/raw/*_prot.mae
for p in /home/users/jwang003/glosa_mods/*
do
	echo $p
	cd $p
	        struct_others=""
		for struct_other in $p/structures/processed/*
		do
			struct_other2=$(basename "${struct_other}")
			#cho $struct2
			struct_others+=" ${struct_other2##*/}"
			job_struct=${struct_other2##*/}
		done
	struct_others=$(echo "$struct_others"|tr ' ' '\n'|tac|tr '\n' ' ')
	echo "$struct_others"
	sbatch -t 8:00:00 -J ${job_struct}_rmsd -o jobs/%j.out -p owners --wrap "python -u /home/users/jwang003/docking/submit_scripts/compute_rmsd_all.py ${struct_others}"
done
