#!/bin/bash

TIMEFORMAT=%R

array=( mpileup_v3.py consensus_v3.py coverage_v3.py) 
for counter in {1..30}
do
	upper_limit=$(($counter*20000))
	
	for command in "${array[@]}"
	do
		echo $command
		data_range="chr3R-${upper_limit}"
		echo $data_range
		start=$(date +%s%N)
		if [ "$command" == 'mpileup' ]
			then
			var1=$(/usr/bin/time -f "%M" -o out.txt  sudo PGPASSWORD=pysam8192 psql -U postgres -h 127.0.0.1 -d pysamdb -c "SELECT public.f_mpileup('chr3R',1,$upper_limit)" )
		elif [ "$command" == 'consensus' ]
			then
			var1=$(/usr/bin/time -f "%M" -o out.txt  sudo PGPASSWORD=pysam8192 psql -U postgres -h 127.0.0.1 -d pysamdb -c "SELECT public.f_consensus('chr3R',1,$upper_limit)" )
		else
			var1=$(/usr/bin/time -f "%M" -o out.txt  sudo PGPASSWORD=pysam8192 psql -U postgres -h 127.0.0.1 -d pysamdb -c "SELECT public.f_coverage('chr3R',1,$upper_limit)" )
		fi
		end=$(date +%s%N)

		var2=`cat out.txt`

		time_range=$((($(date +%s%N) - $start)))
		
		echo $counter
		
		echo "$command,$data_range,$start,$end,$time_range,$var2">> commands_log_PGSQL_10.csv
		
		
	done
	
done
