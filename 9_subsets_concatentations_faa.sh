
#! bin/bash

source /user/work/tt22567/Miniconda3/ENTER/bin/activate   /user/work/tt22567/Miniconda3/ENTER/envs/twist

PWD=/user/work/tt22567/Neuroptera/7_phylogeny/treespace_subsets/subsets_1000_single

cd $PWD

# Define the directory containing the gene files and the folder with the CSV files
gene_file=/user/work/tt22567/Neuroptera/6_treespace/7_treespace/FAA_SUM
csv_folder=./subsets_csv

# Create a directory for the output phylogenetic files
mkdir -p ./subsets_phy

# Loop through each CSV file in the csv_folder
for file in "$csv_folder"/*.csv; do
    
    # Get the current CSV file name, removing the path and suffix
    csv_id=$(basename "$file" .csv)
    
    # Create a directory named subsets_faa with a subdirectory for the current CSV
    subsets_file=./subsets_faa/$csv_id
    mkdir -p $subsets_file
    
    # Read gene names from the second line onwards
    subsets_id=$(tail -n +2 "$file")
    
    # Process each gene
    for i in $subsets_id; do
    
        # Use tr to remove carriage return characters (\r) and sed to remove quotes
        gene_id=$(echo $i | sed -e 's#"##g;s#"##g' | tr -d '\r') &&
    
        # Copy the gene file to the subsets directory
        cp "$gene_file/$gene_id" "$subsets_file"
        
    done &&
    
    # Concatenate the gene files into a single phylogenetic file
    perl catfasta2phyml.pl --concatenate --verbose $subsets_file/*.faa > ./subsets_phy/$csv_id.phy &&
    
    # Remove the temporary subsets directory
    rm -rf $subsets_file
    
done

#Xiumei Lu code