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

# Steps as follows
# 1) Loop the script through a directory/file list of *.spa files - IN PROGRESS
# 2) Copy the first new .csv and rename it `all_spectra.csv` - NOT DONE
# 3) Concatenate all the other files, dropping the wavenumbers - NOT DONE


#### Part 1 - Build the loop


for f in $(find ../data/raw/MIR -name '*.spa')
do
	echo $f
done