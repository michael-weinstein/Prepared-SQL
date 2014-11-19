-- --------------------------
-- Homozygosity identifier --
-- --------------------------

create table {SUBJECT}_homozygosity as select * from
    {SUBJECT}_rare_damage
where
        sampleGenotype NOT IN ('S' , 'R', 'Y', 'M', 'K', 'W') #excludes SNVs with two bases called at the locus
        AND (INSTR(sampleGenotype, '/') = 0 #excludes anything with a slash in the genotype column (as this usually indicates multipe genotypes at the locus)
        OR substring_index(sampleGenotype, '/', 1) = substring_index(sampleGenotype, '/', -1)) #unless the values to the right and left of the slash are the same, because sometimes variant callers put out unexpected stuff
        or (Ref_read_percentVCF < 0.1);  #This looks at the percentage of reads that were reference and can include a variant regardless of what the variant caller called it. My default is 10% or fewer reference reads for possible homozygosity at the locus, but this can be adjusted.
alter table {SUBJECT}_homozygosity ADD INDEX (uniLoc), ADD INDEX (geneList);

/* Usage: This script will generate a new table containing possible homozygous variants in the subject.
As I do not completely trust variant callers to always be right in deciding if an allele is homozygous
or not, I included a line that has it include any position where 90% or more of the reads were
non-reference.  You can change this number if you wish, its selection was arbitrary.
Michael Weinstein, Cohn Lab 2014*/
