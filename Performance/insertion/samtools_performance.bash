#!/bin/bash

TIMEFORMAT=%R

array=(ENCFF469GFN ENCFF074VCI ENCFF194YCE ENCFF324LEM ENCFF920NVC HG002_NA24385_son )

for counter in {1..10}
do
	for file_name in "${array[@]}"
	do
		if [ $file_name != 'HG002_NA24385_son'  ]
		then
			echo 'start $file_name to SAM'
			timing1=$(sudo /usr/bin/time -f "%M" -o out1.txt samtools view -h -o ${file_name}.sam ${file_name}.bam)
			start=$(date "+%Y-%m-%d %H:%M:%S.%N")
			echo ' $file_name to BAM'
			timing1=$(samtools view -bo ${file_name}.bam ${file_name}.sam)
		fi
		if [ $file_name == 'HG002_NA24385_son'  ]
		then
			
			start=$(date "+%Y-%m-%d %H:%M:%S.%N")
		fi
		echo 'start $file_name sorting'
		timing1=$(sudo /usr/bin/time -f "%M" -o out2.txt  samtools sort -o ${file_name}_sorted.bam ${file_name}.bam)
		echo 'start $file_name indexing'
		timing2=$(sudo /usr/bin/time -f "%M" -o out3.txt  samtools index ${file_name}_sorted.bam)
		end=$(date "+%Y-%m-%d %H:%M:%S.%N")
		echo ''
		echo '$file_name ...OK'
		echo ''
		echo $counter
		echo ''
		
		var1=`cat out1.txt`
		var2=`cat out2.txt`
		var3=`cat out3.txt`
		
		memory="${var1} ${var2} ${var3}"
		
		ram_conversion="$(echo "$memory" | cut -d' ' -f1)"
		ram_sorting="$(echo "$memory" | cut -d' ' -f2)"
		ram_indexing="$(echo "$timing" | cut -d' ' -f3)"
		
		find -type f -name '*sorted*' -delete
		
		echo "$counter,$file_name,$start,$end,$ram_conversion,$ram_sorting,$ram_indexing">> sortindex_performance.csv
	done
done
	