####wrote by Xiumei Lu

import csv
import glob

input_dir = './FAA_SUM/'
output_file = '2_output.csv'

# Create an empty list to store the results
results = []

# Loop through all the files in the directory ending with 'bmg.faa'
for filename in glob.glob(input_dir + '*bmg.faa'):

    # Open the file and read the sequences
    with open(filename, 'r') as in_file:
        total_seq = ''
        total_gap = 0
        num_sequences = 0  # Counter for number of fasta sequences
        total_length = 0  # Counter for total length of all sequences

        for i in in_file:
            i=i.strip()   #remove strange strings at tips of sequences
            if i.startswith('>'):
                num_sequences += 1
            else:
                total_seq += i
                total_gap += i.count('-')
                total_length += len(i)

        total_len = len(total_seq)
        total_non_gap = total_len - total_gap
        
        if num_sequences != 0:
            avg_length = total_length / num_sequences
        else:
            avg_length = 0
        
        if total_len != 0:
            gap_ratio = total_gap / total_len
        else:
            gap_ratio = 0
        
        # Append the results to the list
        results.append([filename, num_sequences, total_len, total_gap, total_non_gap, round(gap_ratio, 2), round(avg_length, 2)])


# Open a CSV file in write mode and write the results to the rows
with open(output_file, 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(['Filename', 'Num sequences', 'Total length', 'Total gap', 'Total non-gap', 'Gap ratio', 'Avg length'])
    for row in results:
        writer.writerow(row)

