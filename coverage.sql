
CREATE OR REPLACE FUNCTION public.f_coverage(
	chr_name character varying,
	starting_pos integer,
	ending_pos integer)
    RETURNS TABLE(chromosome_name character varying, start_position integer, end_position integer, numreads integer, covbases bigint, coverage numeric, meandepth numeric, meanbaseq numeric, meanmapq numeric) 
    LANGUAGE 'plpgsql'
    COST 100
    VOLATILE PARALLEL UNSAFE
    ROWS 1000

AS $BODY$

declare counter integer:= (Starting_pos - 50);
declare prev_pos VARCHAR:='';
declare next_pos VARCHAR:='';
declare current_pos VARCHAR:='';
declare is_first_last VARCHAR:='';

declare var_meanmapq NUMERIC(10,2):=0;
declare var_numreads INTEGER:=0;

begin

CREATE TEMPORARY TABLE temp_bas_qual(
	idx INTEGER,
	read_id INTEGER, 
	base_complete VARCHAR,
	base_quality INTEGER,
	read_quality INTEGER
	
);

CREATE TEMPORARY TABLE temp_cigar (
	idx INTEGER,
	pos INTEGER,
	cigar_ext VARCHAR
);

INSERT INTO temp_cigar
select max(b.id),max(b.pos),string_agg(REPEAT(regexp_replace(b.letter,'[0-9]+',''),regexp_replace(b.num,'[A-Z]+','')::INTEGER),'' order by b.idx)
from (
select *,ROW_NUMBER() over() idx
from (
select pos,
regexp_split_to_table(cigar,'(?<=\d+)([A-Z])(?=\d+)') num
,regexp_split_to_table(cigar,'(?<=[A-Z]+)(\d+)(?=[A-Z]+)') letter,
cigar,
id
from Read
where  rname=Chr_name and pos between 	(Starting_pos-50+1) and Ending_pos and cigar <> '*' --Starting_pos and Ending_pos and cigar <> '*'
order by pos) A ) B										
group by b.id
;

--RAISE NOTICE 'counter %, Starting_pos %, Ending_pos %',counter ,Starting_pos,Ending_pos;

while counter < Ending_pos + 1 loop

	insert into temp_bas_qual
	select counter, 
			id read_id,
				CASE
						WHEN substring(cigar_ext from counter+1-Read.pos for 1)='M'
							THEN
								CASE 
	 							WHEN flag=16 THEN LOWER(substring(seq from counter+1-Read.pos +(char_length(replace(substring(cigar_ext from 1 for counter+1-Read.pos),'D',''))-char_length(replace(substring(cigar_ext from 1 for counter+1-Read.pos),'I',''))) for 1)) 
									ELSE UPPER(substring(seq from counter+1-Read.pos +(char_length(replace(substring(cigar_ext from 1 for counter+1-Read.pos),'D',''))-char_length(replace(substring(cigar_ext from 1 for counter+1-Read.pos),'I',''))) for 1))
								END
						WHEN substring(cigar_ext from counter+1-Read.pos for 1)='N' 
							THEN '>'
						ELSE ''
					END base_complete,
				CASE
				WHEN substring(cigar_ext from counter+1-Read.pos for 1)='M'
					THEN ASCII(substring(qual from counter+1-Read.pos +(char_length(replace(substring(cigar_ext from 1 for counter+1-Read.pos),'D',''))-char_length(replace(substring(cigar_ext from 1 for counter+1-Read.pos),'I',''))) for 1))-33
				WHEN substring(cigar_ext from counter+1-Read.pos for 1)='D'
					THEN ASCII(substring(qual from counter+1-Read.pos for 1))-33
				WHEN substring(cigar_ext from counter+1-Read.pos for 1)='I'
					THEN ASCII(substring(qual from counter+1-Read.pos for 1))-33
				ELSE
					0
				END  base_quality,
				Read.mapq quality				
	from Read inner join temp_cigar 
	on Read.id=temp_cigar.idx
	where rname=Chr_name
			and Read.pos between (Starting_pos-50+1) and Ending_pos; --(Starting_pos+1); 

	
	counter := counter + 1;
   	end loop;
	
	--RAISE NOTICE '(Starting_pos-50+1) %, (Starting_pos+1) %, (Ending_pos-Starting_pos-1) %, counter end %',(Starting_pos-50+1), (Starting_pos+1), (Ending_pos-Starting_pos-1),counter;
	
	select count(read_id), ROUND(avg(max_read_quality)::NUMERIC,2) INTO var_numreads,var_meanmapq
	from (
	select read_id,max(read_quality::NUMERIC) max_read_quality
	from temp_bas_qual
	where base_complete <> '' and idx between Starting_pos and Ending_pos
	group by read_id ) a;

	RETURN QUERY
	select Chr_name, Starting_pos,Ending_pos,var_numreads,count(*), ROUND((count(*)/(Ending_pos-Starting_pos+1)::NUMERIC)*100,2),ROUND((sum(pos_count)/(Ending_pos-Starting_pos+1)::NUMERIC),2), ROUND(sum(pos_avg)/sum(pos_count)::NUMERIC,2),var_meanmapq
	from (
	select count(*) pos_count, sum(base_quality::INTEGER) pos_avg
	from temp_bas_qual
	where base_complete <> '' and idx between Starting_pos and Ending_pos
	group by idx
	having count(*)>0) a;
	
	drop table temp_bas_qual;
	drop table temp_cigar;
	

end 
$BODY$;


--example
--select *
--from f_coverage('chr3R',30000000,30000175);

