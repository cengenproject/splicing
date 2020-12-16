#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=spladder-build
#SBATCH -c 17
#SBATCH --mem=10G
#SBATCH --time=3-5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu




bam_dir="/path_to/AB_bss2"
build_dir="/path_to/spladder/spladder_AB_bss2"
ref_gtf="/path_to/references/c_elegans.PRJNA13758.WS274.canonical_geneset.gtf"
tmp_dir="/path_to/tmp/spladder"

mkdir -p $tmp_dir

echo
echo "Running SplAdder on $(date) with setting sg_min_edge_count=3, on the AB_bss2 alignments"
echo


# build list of samples
neuronsList=(ASG AVE AVG AWA AWB PVD VD DD)

for (( i = 0; i < ${#neuronsList[@]}; i++))
do
	sampleArray+=($(ls $bam_dir/${neuronsList[i]}r*.bam ))
done

# Index SAMs if necessary
module load SAMtools
for samp in ${sampleArray[@]}
do
	if [ -f "$samp.bai" ]; then
		echo "$samp already indexed."
	else 
		echo "Indexing $samp"
		samtools index $samp
	fi
done
module purge

function join_with_comma { local IFS=","; echo "$*"; }

comma_bam_list=$(join_with_comma ${sampleArray[@]})

echo "Running command:"
echo
echo "spladder build -o $build_dir \
				-b $comma_bam_list \
				-a $ref_gtf \
				--parallel $SLURM_CPUS_PER_TASK \
				-v \
				--readlen 101 \
				--tmp-dir $tmp_dir \
				--validate-sg \
				--ignore-mismatches"

echo

module load miniconda
source activate spladder 

spladder build -o $build_dir \
				-b $comma_bam_list \
				-a $ref_gtf \
				--parallel $SLURM_CPUS_PER_TASK \
				-v \
				--readlen 101 \
				--tmp-dir $tmp_dir \
				--validate-sg \
				--ignore-mismatches

echo 
echo "Done on $(date)"
echo
