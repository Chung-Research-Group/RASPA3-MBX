#!/bin/bash 

cwdir=$(pwd)

while read struc
do
	while read T
	do
		while read P
		do
			calc_DIR=${cwdir}/calculation/${struc}/${T}K/${P}Pa
			cd $calc_DIR
			if [ -d "output" ]; then
				echo "$calc_DIR pass"
			else
				echo $calc_DIR
				#cp ${cwdir}/raspa3_run.slurm $calc_DIR 
				sbatch -J ${struc}_$P raspa3_run.slurm
			fi
			cd $cwdir
		done < p_range
	done < t_range
done < list

