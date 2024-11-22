import csv
import sys

csv_file = sys.argv[1]
out_file = sys.argv[2]


with open(csv_file, 'r') as f:
    reader = csv.reader(f)
    rows = [row for row in reader]

with open(out_file, 'w') as j:
    for row in rows:
        if row[1] != 'smiles':
            new_row = '\t'.join(row[::-1])
            j.write(new_row + '\n')
