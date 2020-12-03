#!/bin/bash

#***************************************************************#
#   Converts folders of Thermo Omnic *.spa files to one *.csv   #
#																#
#	Written by mark.farrell@csiro.au							#
#***************************************************************#



# Uses a nice ruby script found on https://github.com/ne0dim/spa2csv
# This script assumes the ruby file is called `spa2csv.rb`, and that
# the *.spa files are located in the relative path:
# ../data/raw/MIR
# It assumes there is a directory called ../data/working to put them in
# With time I will try to improve it to automate directory checking etc

## REQUIRES:
# Ruby
# Python
# CSVKit

# Steps as follows
# 1) Loop the script through a directory/file list of *.spa files - DONE
# 2) Concatenate all files, dropping the wavenumbers - DONE
# 3) Clean up

#### Part 0 - Set-up
mkdir -p ../data/working
touch ../data/working/MIRspec.csv

#### Part 1 - Build the loop
for f in $(find ../data/raw/MIR -name '*.spa')
do
	ruby spa2csv.rb $f
done

#### Part 2
csvjoin -c1 $(find ../data/raw/MIR -name '*.csv') > ../data/working/MIRspec.csv

