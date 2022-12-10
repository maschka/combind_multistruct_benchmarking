#!/bin/bash
for p in /home/users/jwang003/glosa_mods/*
do
	echo $p
	cd $p
	for struct in $p/structures/raw/*_prot.mae
	do
		struct2=${struct%_*}
		echo $struct2
		echo "combind structprep ${struct2##*/}"
		sbatch -p rondror -o jobs/%j.out --wrap "combind structprep ${struct2##*/}"
	done
done
