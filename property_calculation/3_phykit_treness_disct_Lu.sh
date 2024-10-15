#! /bin/bash

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/twist

PWD=/user/work/tt22567/Neuroptera/6_treespace/7_treespace

cd $PWD

tree_target=/user/work/tt22567/Neuroptera/6_treespace/7_treespace/MFP_BMGE

rm -rf ./parameters1.csv

if [ ! -f ./parameters1.csv ];then

  touch ./parameters1.csv
  echo -e "tree_id,Treeness,PatristicDist_mean" >> ./parameters1.csv

fi

tree_done=$( awk -F ',' 'NR>1 {print substr($1, index($1,$2))}' parameters1.csv)

for i in $tree_target/*_mft.bmg.faa.treefile; do
  
  id=$( echo $i | sed -e "s#$tree_target\/##g")
  
  Treeness=$(phykit treeness $i)   #1
    
  PatristicDist_mean=$(phykit patristic_distances $i | sed -n 1p | awk '{print $2}')  #8 data，taking the mean

  
  echo -e "$id,$Treeness,$PatristicDist_mean"  >> ./parameters1.csv

done


  