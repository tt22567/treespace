#!/usr/bin/bash

#SBATCH --job-name=1_NS_MAX_1000
#SBATCH --partition=cpu
#SBATCH --account=gely018542
#SBATCH --mem=100GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=7-00:00:0

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/twist

PWD=/user/work/tt22567/Neuroptera/7_phylogeny/treespace_subsets/subsets_1000_single

target_matrix=./subsets_phy/1_NS_MAX_1000.phy   #change here the target

storage_prefix=$(echo $target_matrix | sed -e 's#./subsets_phy/##g;s#.phy##g')        

storage_document=./SimpleModel/$storage_prefix

mkdir -p $storage_document

iqtree -s $target_matrix --score-diff ALL -bb 1000  -m Q.insect+G4  -nt 10 -pre $storage_document/$storage_prefix

#Xiumei Lu code
#take one concatenation alignment as an example