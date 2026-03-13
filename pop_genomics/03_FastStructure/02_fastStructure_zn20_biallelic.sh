#!/bin/bash
#SBATCH --time=0-6:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jessica.l.fletcher@ucdenver.edu
#SBATCH --qos=normal
#SBATCH --partition=amilan

module purge
module load anaconda

conda activate fastStructure_env

mkdir fastStructure_output_zn20_biallelic

for i in {1..10}
do
	python /projects/jfletcher@xsede.org/software/anaconda/envs/fastStructure_env/bin/structure.py -K $i --input=plink_zn20_output --output=fastStructure_output_zn20_biallelic/fastStructure_output_zn20_biallelic
done

python /projects/jfletcher@xsede.org/software/anaconda/envs/fastStructure_env/bin/chooseK.py --input=fastStructure_output_zn20_biallelic/fastStructure_output_zn20_biallelic


python /projects/jfletcher@xsede.org/software/anaconda/envs/fastStructure_env/bin/chooseK.py --input=fastStructure_output_zn20_biallelic/fastStructure_output_zn20_biallelic > fastStructure_output_zn20_biallelic/chooseK_output_zn20_biallelic.txt
