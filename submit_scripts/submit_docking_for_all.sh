#!/bin/bash
for p in /home/users/jwang003/glosa_mods/*
do
  echo $p
	cd $p
	for struct in structures/grids/*
	do
	#struct2=${struct%_*}
		echo $struct
		struct2=${struct##*/}
		echo $struct2
		sbatch -p owners -t 8:00:00 -o jobs/%j.out --wrap "combind dock --grid structures/grids/${struct2}/${struct2}.zip docking ligands/*/*.maegz"
	done
done
