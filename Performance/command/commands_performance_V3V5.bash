#!/bin/bash

TIMEFORMAT=%R

#touch commands_log.csv
array=( mpileup_v3.py consensus_v3.py coverage_v3.py) 
for counter in {1..10}
do
	for counter in {1..30}
	do
		upper_limit=$(($counter*20000))
		
		for command in "${array[@]}"
		do
			start=$(date +%s%N)
	
			var1=$(/usr/bin/time -f "%M" -o out1.txt  python3 $command -r chr3R:1-$upper_limit -p )
	
			end=$(date +%s%N)
	
			var2=`cat out1.txt`
			timing="${var1} ${var2}"
	
			time_range=$((($(date +%s%N) - $start)))
			data_range="chr3R-${upper_limit}"
			db_min="$(echo "$timing" | cut -d' ' -f1)"
			overall_min="$(echo "$timing" | cut -d' ' -f2)"
			ram_command="$(echo "$timing" | cut -d' ' -f3)"
			
			echo $counter
			
			echo "$command,$data_range,$start,$end,$db_min,$overall_min,$ram_command">> coverage_final_V3V5.csv
			
			
		done
		
	done
done
