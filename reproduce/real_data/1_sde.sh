#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --job-name=1_sde
#SBATCH --partition=main,mulan
#SBATCH --mem=50G
#SBATCH --array=1-626
#SBATCH --output=/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/out/out_%a.out
#SBATCH --error=/net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/out/out_%a.err

bash

let k=0
awk 'NR > 1 {print $2}' /net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/trait.txt | while read value; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
Rscript 1_sde.R $value 
if [ -e "$value/GWAS.RDat" ]; then
rm /net/mulan/disk2/luliuu/project4/realdata/Gtex/GWASsex5e_8/out/out_${SLURM_ARRAY_TASK_ID}.*
fi
fi
done
