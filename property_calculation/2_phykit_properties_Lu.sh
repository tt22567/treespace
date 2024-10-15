#! /bin/bash

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/twist

PWD=/user/work/tt22567/Neuroptera/6_treespace/7_treespace

cd $PWD

gene_target=/user/work/tt22567/Neuroptera/6_treespace/7_treespace/FAA_SUM
tree_target=/user/work/tt22567/Neuroptera/6_treespace/7_treespace/TREEFILE_MFP_SUM
out_group=/user/work/tt22567/Neuroptera/6_treespace/7_treespace/outgroup.csv

gene_id=$( echo $tree_target/*mft.bmg.faa.treefile | sed -e "s#$tree_target\/##g;s#.treefile##g") #no empty genes


#rm -rf ./parameters.csv

if [ ! -f ./parameters.csv ];then

  touch ./parameters.csv
  echo -e "gene_id,GCcontent(%),PairwiseIden_mean,ParsiInfoSite_Num,ParsiInfoSite_ratio,RelaCompVar,VariSite_Num,VariSite_ratio,tree_id,bipss_mean,dvMolClock,EvoRat,LonBranScore_mean,TotalTreeLen,Treeness,PatristicDist_mean,Saturation,TreenessdeRcv,date" >> ./parameters.csv

fi

gene_done=$( awk -F ',' 'NR>1 {print substr($1, index($1,$2))}' parameters.csv)

for i in $gene_id; do
  
  if [[ ! $gene_done =~ $i ]]; then
  
  tree_id=${i}.treefile

  #for genes
  GCcontentPerc=$(phykit gc_content $gene_target/$i)
  
  PairwiseIden_mean=$(phykit pairwise_identity  $gene_target/$i | sed -n 1p | awk '{print $2}' )  #8 data, taking the mean
  
  ParsiInfoSite_Num=$(phykit parsimony_informative_sites $gene_target/$i | sed " " | awk '{print $1}' ) #3 data, taking the Number
  ParsiInfoSite_ratio=$(phykit parsimony_informative_sites $gene_target/$i | sed " " | awk '{print $3}' ) #3 data, taking the ratio
  
  RelaCompVar=$(phykit relative_composition_variability $gene_target/$i )  
  
  VariSite_Num=$(phykit variable_sites $gene_target/$i | sed " " | awk '{print $1}' ) #3 data, taking the Number
  VariSite_ratio=$(phykit variable_sites $gene_target/$i | sed " " | awk '{print $1}' ) #3 data, taking the ratio
  
  #for tree
  bipss_mean=$(phykit bipartition_support_stats $tree_target/$tree_id | sed -n 1p | awk '{print $2}' ) #8 data, taking the mean
  
  dvMolClock=$(phykit degree_of_violation_of_a_molecular_clock -t $tree_target/$tree_id -r $out_group )
  
  EvoRat=$(phykit evolutionary_rate $tree_target/$tree_id) 
  
  LonBranScore_mean=$(phykit long_branch_score $tree_target/$tree_id | sed -n 1p | awk '{print $2}' ) #8 data, taking the mean value
  
  TotalTreeLen=$(phykit total_tree_length $tree_target/$tree_id)
  
  Treeness=$(phykit treeness $tree_target/$tree_id)   #
  
  PatristicDist_mean= $(phykit patristic_distances $tree_target/$tree_id | sed -n 1p | awk '{print $2}')  #8 data, taking the mean value

  #for gene and tree
  Saturation=$(phykit saturation -a $gene_target/$i  -t $tree_target/$tree_id)
  
  TreenessdeRcv=$(phykit treeness_over_rcv -a $gene_target/$i  -t $tree_target/$tree_id |  sed " " | awk '{print $3}' ) #3 data, taking the Treeness/RCV
  
  
  echo -e "$i,$GCcontentPerc,$PairwiseIden_mean,$ParsiInfoSite_Num,$ParsiInfoSite_ratio,$RelaCompVar,$VariSite_Num,$VariSite_ratio,$tree_id,$bipss_mean,$dvMolClock,$EvoRat,$LonBranScore_mean,$TotalTreeLen,$Treeness,$PatristicDist_mean,$Saturation,$TreenessdeRcv,$(date)"  >> ./parameters.csv

  fi 

done

#GC content is not suitable for protein sequences