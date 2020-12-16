#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=spladder-test
#SBATCH -c 15
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL


module load miniconda
source activate spladder 

bam_dir="/path_to/AB_bss2"
build_dir="/path_to/spladder/spladder_AB_bss2"
ref_gtf="/path_to/references/c_elegans.PRJNA13758.WS274.canonical_geneset.gtf"
tmp_dir="/path_to/tmp/spladder"


echo
echo "Running SplAdder test for all pairs of neurons (8 types) on $(date) with setting sg_min_edge_count=3, on the AB_bss2 alignments"
echo



neuronsList=(ASG AVE AVG AWA AWB PVD VD DD)





for (( i = 0; i < ${#neuronsList[@]}; i++))
do
	for (( j = $i+1; j < ${#neuronsList[@]}; j++))
	do
		echo "##################################################################"
		echo "i: "$i", j: "$j
		neurA=${neuronsList[$i]}
		neurB=${neuronsList[$j]}
		
		bamsA=$(ls $bam_dir/${neuronsList[i]}r*.bam | tr '\n' ',' | sed 's/,$//')
		bamsB=$(ls $bam_dir/${neuronsList[j]}r*.bam | tr '\n' ',' | sed 's/,$//')
		
		echo "Testing $neurA vs $neurB"
		echo
		echo "Full command:
		spladder test --conditionA $bamsA \
--conditionB $bamsB \
--outdir $build_dir \
--parallel $SLURM_CPUS_PER_TASK \
--readlen 101 \
--validate-sg \
--labelA $neurA --labelB $neurB \
--diagnose-plots"
			  
		spladder test --conditionA $bamsA \
					  --conditionB $bamsB \
					  --outdir $build_dir \
					  --parallel $SLURM_CPUS_PER_TASK \
					  --readlen 101 \
					  --validate-sg \
					  --labelA $neurA --labelB $neurB \
					  --diagnose-plots

		echo
		echo
	done
done



echo 
echo "Done on $(date)"


