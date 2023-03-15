#!/usr/bin/env python3

'''
This script converts an asc file to a list of pixel coordinates for presences.
Assumes the input is a thresholded binary map; missing data values don't matter as long as they aren't 1.

Example loop:

for f in combined_used_ranges/*.asc; do
i="$(echo ${f} | sed 's/\.asc//g' | sed 's/.*\///g')"
echo ${i}
./asc_to_PAM.py ${f} ${i}.pam
done

When this is done, finish by running:
cat *.pam > pam_all.pam
sort pam_all.pam > temp && mv temp pam_all.pam
mkdir pam_files; mv *.pam pam_files/; mv pam_files/pam_all.pam ./

Then run combine_PAM.py

'''

import csv
import numpy 
import sys 
import argparse # Parse arguments
import subprocess
import re
import random
import scipy.spatial
import itertools

parser = argparse.ArgumentParser(description='Script to delete unnecessary fields in locality documets that are tab-delimited and follow DarwinCore header terms.')
parser.add_argument('input_file', action='store', help='Name of the asc input file.')
parser.add_argument('output_file', action='store', help='Name of the desired output file.')
args = parser.parse_args()

infile = args.input_file # Input asc
outfile = args.output_file # Output 

#for file in infile:
matrix = []
with open(infile, 'r') as datafile:
	reader=csv.reader(datafile, delimiter=' ')
	for _ in range(6):
		next(reader) # Skip first 6 rows
	for r in reader:
		r = list(filter(None, r)) # Remove empty cells -- needed for leading spaces
		matrix.append(r)		

cellsize = float(subprocess.check_output("grep cellsize {} | sed 's/.* //g'".format(infile), shell = True))

xcorner = float(subprocess.check_output("gdalinfo {} | grep 'Upper Left' | sed 's/,.*//g' | sed 's/.*(//g'".format(infile), shell = True))
xcorner = xcorner + 180 # Make longitude -180 become 0
ycorner = float(subprocess.check_output("gdalinfo {} | grep 'Upper Left' | sed 's/).*//g' | sed 's/.*, *//g'".format(infile), shell = True))
ycorner = 90 - ycorner # Make latitude 90 become zero; equator would be 90
print(ycorner)


xcorner_integer = int(round(xcorner/cellsize))
ycorner_integer = int(round(ycorner/cellsize))

print(infile)
species = re.sub('.*\/', '', infile)
print(species)
species = re.sub('\..*', '', species)
print(species)

print(", ".join(map(str, [cellsize, xcorner, ycorner]))) # Have to convert list of floats to list of strings for join method
print(", ".join(map(str, [xcorner_integer, ycorner_integer])))

# Use these two variables to track our coordinates in integers
x_coord = xcorner_integer
y_coord = ycorner_integer

result = []
# Loop over the dataset
for r in matrix:
	for i in r:
		if i == '1':
			result.append(["_".join([str(x_coord), str(y_coord)]), species])
		x_coord += 1
	y_coord += 1 # Row number corresponds to y dimension numbered in pixels
	x_coord = xcorner_integer # We have finished this row, so we reset the X coordinate to read through the next one

#print(result)
	
with open(outfile, 'w+') as writefile:
	writer = csv.writer(writefile, delimiter='\t')
	for r in result:
		writer.writerows([r]) # The syntax on this line was problematic in terms of iteration behavior


