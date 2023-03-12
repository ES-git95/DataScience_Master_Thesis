#!/bin/bash

TIMEFORMAT=%R

#touch commands_log.csv
for iteration in {1..10}
do
	array=( mpileup consensus coverage) 
	for counter in {1..30}
	do
		upper_limit=$(($counter*20000))

		for command in "${array[@]}"
		do
			start=$(date +%s%N)

			if [ "$command" == "mpileup" ]
				then
				var1=$(/usr/bin/time -f "%M" -o out1.txt  samtools mpileup ENCFF469GFN_sorted.bam -r chr3R:1-$upper_limit )
			elif [ "$command" == "consensus" ]
				then
				var1=$(/usr/bin/time -f "%M" -o out1.txt  samtools consensus ENCFF469GFN_sorted.bam -r chr3R:1-$upper_limit --mode simple --call-frac 0.55 )
			else 
				var1=$(/usr/bin/time -f "%M" -o out1.txt  samtools coverage ENCFF469GFN_sorted.bam -r chr3R:1-$upper_limit )
			fi
			echo ""	

			end=$(date +%s%N)

			var2=`cat out1.txt`
			timing="${var2}"

			time_range=$((($(date +%s%N) - $start)))
			data_range="chr3R-${upper_limit}"
			ram_command="$(echo "$timing" | cut -d' ' -f1)"

			echo $counter

			echo "$command,$data_range,$start,$end,$time_range,$ram_command">> commands_log_samtools_final.csv


		done

	done
done