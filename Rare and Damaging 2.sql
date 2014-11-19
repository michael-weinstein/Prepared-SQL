-- ------------------------------------------------------------------------------------------------------------
-- This creates a new table of only rare mutations with a decent chance of being damaging to the gene product--
-- ------------------------------------------------------------------------------------------------------------

 
create table {SUBJECT}_rare_damage as select *
from
    {SUBJECT}

left join {SUBJECT}VCF on uniLoc = uniLocVCF
where
    (AfricanHapMapFreq < 1.0)
	and (AsianHapMapFreq < 1.0)
	and (EuropeanHapMapFreq < 1.0)
	and (dbSNP_5_percent_in_all != 'TRUE')
	and (dbSNP_5_percent_in_any != 'TRUE')
	and (MAFinESP < 0.01) #selects only variants with frequency < 1% in each database (or <5% in every dbSNP population
	and (functionGVS NOT IN ('intron' , 'synonymous',
        '5-prime-UTR',
        '3-prime-UTR',
        'non-coding-exon',
        'upstream-gene',
        'downstream-gene',
        'intergenic')) group by uniLoc;   #deselects any variants annotated with the listed consequences
alter table {SUBJECT}_rare_damage add primary key (lineNumber, uniLoc, chromosome, position, geneList);


/*Usage of this code: This will generate a new table containing the rare mutations which are likely to affect
the function of the gene product.  Please replace {SUBJECT} with the actual name of your table, otherwise this will
(obviously) fail.  This snippet of code was written assuming that you ran both dbSNPerator2.5.pl and 
seattleSequel1.1.pl on your output.  If that is not the case, you will need to change some column names and
probably eliminate a few of them.

Michael Weinstein, Cohn Lab 2014*/