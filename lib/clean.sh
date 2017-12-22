#!/bin/bash
## Author: A.Hernandez
## version v2.0

if [ $# -eq 0  ]; then
        echo -e "\nExecute cleaning\u"
        echo "Usage: clean.sh <config.file>"
        exit
fi

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

#VARIABLES
dir=$1

if [ $mapping == "YES" ]; then

	if [ $duplicate_filter == "YES" ]; then
		find $dir/Alignment -name *align* | xargs -I % echo "rm %" >> $dir/logs/cleanning.txt
	else
		find $dir/Alignment -name *align.b* | xargs -I % echo "rm %" >> $dir/logs/cleanning.txt
	fi
fi

if [ $cfsan == "YES" ]; then
	find $dir/CFSAN/samples -name reads.all.pileup | xargs -I % echo "rm %" >> $dir/logs/cleanning.txt
	find $dir/CFSAN/samples -name reads.sam | xargs -I % echo "rm %" >> $dir/logs/cleanning.txt
	find $dir/CFSAN/samples -name reads.sorted.bam | xargs -I % echo "rm %" >> $dir/logs/cleanning.txt
	find $dir/CFSAN/samples -name reads.unsorted.bam | xargs -I % echo "rm %" >> $dir/logs/cleanning.txt
fi

if [ $srst2 == "YES" ]; then
	find $dir/srst2 -name *.pileup | xargs -I % echo "rm %" >> $dir/logs/cleanning.txt
	find $dir/srst2 -name *.sorted.bam | xargs -I % echo "rm %" >> $dir/logs/cleanning.txt
fi

bash $dir/logs/cleanning.txt
