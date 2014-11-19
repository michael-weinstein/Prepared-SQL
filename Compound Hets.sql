-- ---------------------------------------------------
-- Find genes with possible compound heterozygosity --
-- ---------------------------------------------------

create table {SUBJECT}_hit_count as select distinct uniLoc, geneList, count(*) hits  #creates a new table containing the listed fields plus a new count field that will be aliased as "hits"
from
    {SUBJECT}_rare_damage
group by geneList;  #tells SQL to group and count by genelist (return entries per genelist value)

create table {SUBJECT}_compound_het as select * from
    {SUBJECT}_rare_damage
where
    geneList in (select distinct
            geneList
        from
            {SUBJECT}_hit_count
        where
            hits > 1) #creates a list of genes with more than one hit from the hit count table and adds them to the new list if the count is more than one
order by geneList asc;  

/*drop table {SUBJECT}_hit_count;*/ #if you really, really need to

/*Usage: This script will first generate a table counting the number of rare, potentially damaging variants (hits) 
there are in each gene.  This script will then use the data in that table to generate a new table listing variants
possibly contributing to a condition of compound heterozygosity for the listed genes using the data in the
previously-generated table of hit counts.  The commented final line will delete the table of hit counts after 
generating the compound heterozygosity table.  I do not recommend this, as it unnecessarily deletes information, but
it should not change your final results.
This script assumes you have already run the script to isolate the rare and potentially damaging
variants in their own table called {SUBJECT}_rare_damage.  If you have not, this will return several genes as
potential compound heterozygosity based on common variants, synonymous SNPs, and the like (and then only if you
set it to run on your unfiltered original table).  

Michael Weinstein, Cohn Lab 2014*/