#!/bin/bash
## Author: A.Hernandez
## Usage: summary_srst2.sh ...


# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x


# VARIABLES

dir=$1
output_dir=$2


echo -e "Finding resistence results"
find $dir -name "*__genes__ARGannot.r1__results.txt" |xargs -I % echo "ln -s % ."| bash

echo -e "Finding mlst results"
find $dir -name  "*mlst__mlst__*__results.txt" |xargs -I % echo "ln -s % ."| bash

echo -e "Running summary_srst2.sh"
srst2 --prev_output *.txt --output $output_dir/srst2

echo -e "Delete ln"
rm *results.txt

echo -e "summary_srst2.sh finished"
