-- FUNCTION: public.f_consensus(character varying, integer, integer)

-- DROP FUNCTION IF EXISTS public.f_consensus(character varying, integer, integer);

CREATE OR REPLACE FUNCTION public.f_consensus(
	chr_name character varying,
	starting_pos integer,
	ending_pos integer)
    RETURNS void
    LANGUAGE 'plpgsql'
    COST 100
    VOLATILE PARALLEL UNSAFE
AS $BODY$

declare counter integer:=Starting_pos;
declare count_A INTEGER:=0;
declare count_G INTEGER:=0;
declare count_T INTEGER:=0;
declare count_C INTEGER:=0;

declare num_rows integer:=0;
declare read_length integer:=0;

begin

CREATE TEMPORARY TABLE temp_bas_qual(
	idx INTEGER,
	read_id INTEGER,
	base_complete char(1)
	
);

CREATE TEMPORARY TABLE consensus_bases(
	base VARCHAR
	
);

CREATE TEMPORARY TABLE consensus_seq(
	rname VARCHAR,
	start_pos INTEGER,
	end_pos INTEGER,
	consensus VARCHAR
	
);

select length(seq) into read_length
from Read
limit 1;

while counter < Ending_pos+1 loop

	insert into temp_bas_qual
	select counter, 
			Read.id read_id,
				case substring(Read.seq from counter+1-Read.pos for 1)
				when '' then ' '
				else substring(Read.seq from counter+1-Read.pos for 1)
				end
				base_complete				
	from Read 
	where Read.rname=Chr_name
			and Read.pos between Starting_pos-read_length-1 and Ending_pos; 		
	
	counter := counter + 1;
	end loop;
	
insert into consensus_bases(base)
SELECT substring(string_agg(base_complete,'' ORDER BY total DESC),1,1)
FROM (
select idx, base_complete,count(*) TOTAL
from temp_bas_qual
group by idx,base_complete
order by idx 
	) TEMP
group by idx;

select count(*) into num_rows
from consensus_bases;

counter:=0;

while counter < num_rows loop
	insert into consensus_seq(rname,start_pos,end_pos,consensus)
	select Chr_name, 
		Starting_pos + counter,
		Starting_pos + counter + read_length-1,
		substring(string_agg(case when base = '' then 'N' ELSE base END,'') from 1 for (Ending_pos-counter))
	from (
		select base
		from consensus_bases 
		limit read_length OFFSET counter) a; 		

counter := counter + read_length;
end loop;

--RETURN QUERY 
--select *
perform *
from consensus_seq;

drop table consensus_bases;
drop table consensus_seq;
drop table temp_bas_qual;
	
	
end 
$BODY$;

ALTER FUNCTION public.f_consensus(character varying, integer, integer)
    OWNER TO postgres;
