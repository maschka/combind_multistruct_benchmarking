#!/bin/bash
#for struct in ../data/*/structures/raw/*_prot.mae
for p in /home/users/jwang003/glosa_mods/*
do
	echo $p
	cd $p
	# sbatch -p owners --wrap "combind ligprep raw/additional*.smi ligands"
  sbatch -p owners -o jobs/%j.out --wrap "combind ligprep raw/*native.smi ligands"
done
