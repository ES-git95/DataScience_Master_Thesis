#!/bin/bash

TIMEFORMAT=%R

array=( ENCFF469GFN.bam ENCFF272TLP.bam ENCFF920NVC.bam ENCFF324LEM.bam HG002_NA24385_son.bam) 


version=1


for file in "${array[@]}"
do
	echo $file
	echo "dropper started..."
	if [ "$version" == 1 ]
		then
		echo "skipped"
		#python3 DBdeleater_pg4.py
	else
		rm PysamDB.db
	fi
	echo "dropper OK"
	echo ""
	echo "creator started..."
	if [ "$version" == 5 ]
		then
		#echo "skipped"
		python3 1_DBcreator_v5.py
	elif [ "$version" == 3 ]
		then
		python3 1_DBcreator_v3.py
	else 
		python3 DBcreator_v3_pg.py
	fi
	echo "creator OK"
	echo ""
	echo "insertion started..."
	echo ""
	if [ "$version" == 5 ]
		then
		start=$(date +%s%N)
		var1=$(/usr/bin/time -f "%M" -o out1.txt  python3 2_DBinsertion_v5.py -f $file -p )
		#timing="$(sudo python3 2_DBinsertion_v5.py -f $file -p)"
		end=$(date +%s%N)
	elif [ "$version" == 3 ]
		then
		start=$(date +%s%N)
		var1=$(/usr/bin/time -f "%M" -o out1.txt  python3 2_DBinsertion_v3.py -f $file -p )
		#timing="$(python3 2_DBinsertion_v3.py -f $file -p)"
		end=$(date +%s%N)
	else 
		start=$(date +%s%N)
		var1=$(/usr/bin/time -f "%M" -o out1.txt  python3 DBinsertion_v4_pg.py -f $file -p )
		#timing="$(python3 DBinsertion_v4_pg.py -f $file -p)"
		end=$(date +%s%N)
	fi
	echo ""
	echo "insertion OK"
	if [ "$version" == 5 ]
		then
		var2=$(/usr/bin/time -f "%M" -o out2.txt  python3 3_DBindex_creator_v5.py )
		#python3 3_DBindex_creator_v5.py
		after_index=$(date +%s%N)
		echo "index creation executed"
	elif [ "$version" == 3 ]
		then
		var2=$(/usr/bin/time -f "%M" -o out2.txt  python3 3_DBindex_creator_v3.py )
		#python3 3_DBindex_creator_v3.py
		after_index=$(date +%s%N)
		echo "index creation executed"
	else 
		var2=$(/usr/bin/time -f "%M" -o out2.txt  python3 DBindex_creator.py )
		#python3 DBindex_creator.py
		after_index=$(date +%s%N)
		echo "index creation executed"
		echo ""
		read -p "please drop tables on database and press a key to continue..."
		echo ""	
	fi
	
	var3=`cat out1.txt`
	var4=`cat out2.txt`
	timing="${var1} ${var3} ${var4}"
	
	time_range=$(($end - $start))
	time_after_index=$(($after_index - $start))
	db_min="$(echo "$timing" | cut -d' ' -f1)"
	overall_min="$(echo "$timing" | cut -d' ' -f2)"
	ram_insertion="$(echo "$timing" | cut -d' ' -f3)"
	ram_index="$(echo "$timing" | cut -d' ' -f4)"
		
	echo $overall_min
	echo ""
	echo ""
	echo "$version,$file,$start,$end,$after_index,$time_range,$time_after_index,$db_min,$overall_min,$ram_insertion,$ram_index,$timing">> insertions_log_PGSQL.csv
	
done
#done
