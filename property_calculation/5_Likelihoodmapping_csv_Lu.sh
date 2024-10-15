#! bin/bash
#SBATCH --job-name=LHM
#SBATCH --account=gely018542
#SBATCH --partition=cpu
#SBATCH --mem=80GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=10
#SBATCH --time=14-00:00:00
#SBATCH --output=./bin1_1_genes_%A.out

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/twist

PWD=/user/work/tt22567/Neuroptera/6_treespace/7_treespace

cd $PWD

storage_file=./MFP_BMGE

touch ./likelihood.csv

echo 'gene_id,fully resolve, partly resolved, unresolved,date' >> ./likelihood.csv

for i in $storage_file/*.lmap.eps;do

  id=$(echo $i | sed "s#$storage_file\/##g;s#.lmap.eps##g")
  od=${id}.iqtree

  value1=$(cat $storage_file/$od | grep 'Number of fully resolved  quartets (regions 1+2+3):' | awk -F "=" '{print $2}' | awk -F "%" '{print $1}')
  value2=$(cat $storage_file/$od | grep 'Number of partly resolved quartets (regions 4+5+6):' | awk -F "=" '{print $2}' | awk -F "%" '{print $1}')
  value3=$(cat $storage_file/$od | grep 'Number of unresolved      quartets (region 7)     : ' | awk -F "=" '{print $2}' | awk -F "%" '{print $1}')
  
  echo -e "$od,$value1,$value2,$value3,${date}" >> ./likelihood.csv
done





  