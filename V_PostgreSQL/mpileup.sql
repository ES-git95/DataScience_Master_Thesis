-- FUNCTION: public.f_mpileup(character varying, integer, integer)

-- DROP FUNCTION IF EXISTS public.f_mpileup(character varying, integer, integer);

CREATE OR REPLACE FUNCTION public.f_mpileup(
	chr_name character varying,
	starting_pos integer,
	ending_pos integer)
    RETURNS void
    LANGUAGE 'plpgsql'
    COST 100
    VOLATILE PARALLEL UNSAFE
AS $BODY$

declare counter integer:=	Starting_pos;			
declare prev_pos VARCHAR:='';
declare next_pos VARCHAR:='';
declare current_pos VARCHAR:='';
declare is_first_last VARCHAR:='';
declare read_length integer:=0;

begin

CREATE TEMPORARY TABLE temp_bas_qual(
	idx INTEGER,
	bases VARCHAR, 
	qualities VARCHAR,
	pos INTEGER
);

CREATE TEMPORARY TABLE temp_cigar (
	idx INTEGER,
	pos INTEGER,
	cigar_ext VARCHAR
);

select length(seq) into read_length
from Read
limit 1;

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
where  rname=Chr_name and pos between 	(Starting_pos-read_length+1) and Ending_pos and cigar <> '*' 
order by pos) A ) B										
group by b.id
;

while counter < Ending_pos +1 loop

	insert into temp_bas_qual
	select counter, 
				CASE 
				WHEN substring(cigar_ext from counter+1-Read.pos for 1) ='D'
					THEN '*'
				ELSE 
					CONCAT(CASE
							WHEN counter=Read.pos
								THEN concat('^',chr(Read.mapq + 33))
							ELSE ''
					END,
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
					END , CASE
							WHEN counter=Read.pos+length(cigar_ext)-1
								THEN '$'
							ELSE ''
					END,CASE 
						WHEN substring(cigar_ext from counter+2-Read.pos for 1)='D' and substring(cigar_ext from counter+1-Read.pos for 1)<>'D'
							THEN CASE
						   		WHEN flag=16 then CONCAT('-', substring(cigar from POSITION('D' in cigar)-1 for 1),'n')
						   		WHEN flag<>16 then CONCAT('-', substring(cigar from POSITION('D' in cigar)-1 for 1),'N')
						   		END
						WHEN substring(cigar_ext from counter+2-Read.pos for 1)='I' and substring(cigar_ext from counter+1-Read.pos for 1)<>'I'
							THEN 
						   		CASE
						   		WHEN flag=16 THEN CONCAT('+', substring(cigar from POSITION('I' in cigar)-1 for 1),REPEAT(LOWER(substring(seq from counter+2-Read.pos for 1)),substring(cigar from POSITION('I' in cigar)-1 for 1)::INTEGER))
						   		ELSE CONCAT('+', substring(cigar from POSITION('I' in cigar)-1 for 1),REPEAT(UPPER(substring(seq from counter+2-Read.pos for 1)),substring(cigar from POSITION('I' in cigar)-1 for 1)::INTEGER))
						   		END
						   END,
						   CASE
						   WHEN substring(cigar_ext from counter+1-Read.pos for 1)='I' and flag=16
						   	THEN lower(substring(seq from counter+1-Read.pos +(char_length(replace(substring(cigar_ext from 1 for counter+1-Read.pos),'D',''))-char_length(substring(cigar_ext from 1 for counter+1-Read.pos))) for 1))
						   WHEN substring(cigar_ext from counter+1-Read.pos for 1)='I' and flag<>16
						   	THEN upper(substring(seq from counter+1-Read.pos +(char_length(replace(substring(cigar_ext from 1 for counter+1-Read.pos),'D',''))-char_length(substring(cigar_ext from 1 for counter+1-Read.pos))) for 1))
						   ELSE ''
						   END
					)
				END base,
				CASE
				WHEN substring(cigar_ext from counter+1-Read.pos for 1)='M'
					THEN substring(qual from counter+1-Read.pos +(char_length(replace(substring(cigar_ext from 1 for counter+1-Read.pos),'D',''))-char_length(replace(substring(cigar_ext from 1 for counter+1-Read.pos),'I',''))) for 1)
				WHEN substring(cigar_ext from counter+1-Read.pos for 1)='D'
					THEN substring(qual from counter+1-Read.pos for 1)
				WHEN substring(cigar_ext from counter+1-Read.pos for 1)='I'
					THEN substring(qual from counter+1-Read.pos for 1)
				ELSE
					''
				END 
				quality,
				Read.pos
	from Read inner join temp_cigar 
	on Read.id=temp_cigar.idx
	where rname=Chr_name
			and Read.pos between (Starting_pos-read_length+1) and LEAST(counter+read_length,Ending_pos) 
	order by Read.pos;
	
	counter := counter + 1;
   	end loop;
	

	perform Chr_name,idx,'N',length(replace(string_agg(qualities,''),'-','')),string_agg(bases,''),replace(string_agg(qualities,''),'-','') 
	from temp_bas_qual
	where qualities <> ''
	group by idx
	order by idx;
	
	drop table temp_bas_qual;
	drop table temp_cigar;

end 
$BODY$;

ALTER FUNCTION public.f_mpileup(character varying, integer, integer)
    OWNER TO postgres;
