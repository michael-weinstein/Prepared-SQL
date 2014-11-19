select round(least(depth_proband, depth_mom, depth_dad), -1) bin,
       count(*) n
from {FAMILY}_genotypes
where  #this is a somewhat complex query explained in detail in the comment below
(substring_index(proband_geno, '/', 1) != substring_index(proband_geno, '/', -1)) and  #tests to exclude loci that are homozygous in proband
((substring_index(proband_geno, '/', 1) != '0' and substring_index(proband_geno, '/', 1) not in #excludes an allele that is a 0 (reference read) from further analysis
  (substring_index(mom_geno, '/', 1), substring_index(mom_geno, '/', -1),  #this line and the next see if the observed allele is seen in mom or dad
   substring_index(dad_geno, '/', 1), substring_index(dad_geno, '/', -1))
) or  
(substring_index(proband_geno, '/', -1) != '0' and substring_index(proband_geno, '/', -1) not in  #does the same as above for the proband's other allele
  (substring_index(mom_geno, '/', 1), substring_index(mom_geno, '/', -1),
   substring_index(dad_geno, '/', 1), substring_index(dad_geno, '/', -1))
))
 group by bin;