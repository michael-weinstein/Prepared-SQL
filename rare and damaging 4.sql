-- -------------------------------------------------------
-- Firebolt Rare and potentially damaging variant finder--
-- -------------------------------------------------------
#alter table {SUBJECT}VCF add column uniqueIDVCF varchar(255);
#update {SUBJECT}VCF set uniqueIDVCF = concat(unilocVCF,altVCF);
create table {SUBJECT}_rare_damage_temp as select *
from
    {SUBJECT}
where
    (AfricanHapMapFreq < 1.0)
	and (AsianHapMapFreq < 1.0)
	and (EuropeanHapMapFreq < 1.0)
	and (dbSNP_5_percent_in_all != 'TRUE')
	and (dbSNP_5_percent_in_any != 'TRUE')
	and (F_rarest_allele < 0.01)  #selects only variants where no population in ExAC has a fequency over 1%
	and (MAFinESP < 0.01) #selects only variants with frequency < 1% in each database (or <5% in every dbSNP population
	and (functionGVS NOT IN ('intron' , 'synonymous',
        '5-prime-UTR',
        '3-prime-UTR',
        'non-coding-exon',
        'upstream-gene',
        'downstream-gene',
        'intergenic')) group by uniLoc;   #deselects any variants annotated with the listed consequences
delete from {SUBJECT}VCF where instr(uniLocVCF, 'GL') != 0;  #removes from the VCF table any lines that do not correspond to a chromosome
#ALTER TABLE {SUBJECT}VCF MODIFY uniLocVCF varchar(15) NOT NULL UNIQUE;  #declares to mySQL that the VCF locus column will have an entry for every line and that entry will be unique
#Alter table {SUBJECT}_rare_damage_temp add CONSTRAINT `uniqueIDVCF_108_0`  #adds a foreign key telling the program that the uniLoc (locus) value will link back to a locus value in the VCF table
#    FOREIGN KEY (`uniqueIDVCF`)
#    REFERENCES `bo`.`{SUBJECT}VCF` (`uniLocVCF`);
create table {SUBJECT}_rare_damage as select * from {SUBJECT}_rare_damage_temp left join {SUBJECT}VCF on uniLoc = uniLocVCF;  #creates a new table by combining the data from the temporary table generated at the beginning of this script with the locus-matched data from the VCF table of read data
drop table {SUBJECT}_rare_damage_temp;  #now that we've created the final table with all the data, we will delete the temporary table with only half the data
ALTER TABLE {SUBJECT}_rare_damage ADD INDEX (lineNumber), ADD INDEX (uniLoc), ADD INDEX (geneList);
/*alter table {SUBJECT}_rare_only add index (uniLoc);  #adds indexing to some lines we may later want to search on
alter table {SUBJECT}_rare_only add index (geneList);*/