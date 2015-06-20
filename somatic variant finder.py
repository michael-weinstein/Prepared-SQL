class SomaticChanges(object):
    
    def __init__(self, proband, parent1 = "", group = ""):
        if not (proband and parent1) or not type(proband) is str:
            try:
                assert hasattr(proband, '__iter__') and len(proband) == 2
                self.group = parent1
                self.parent1 = proband[1]
                self.proband = proband[0]
            except:
                raise RuntimeError("Finding changed variants requires passing proband and original as either arguments or a tupple of two identifiers")
        else:
            self.proband = proband
            self.group = group
            self.parent1 = parent1
            
        
    def buildTables(self):
        query = ""        
        query += "DELETE FROM " + self.proband + "VCF where instr(uniLocVCF, 'GL') != 0;\n"
        query += "DELETE FROM " + self.parent1 + "VCF where instr(uniLocVCF, 'GL') != 0;\n"
        
        query += "CREATE TABLE " + self.proband + self.group + "_genotemp1 as select uniLoc, varMD5, geneList, GenotypeVCF proband_geno, \
\nfunctionGVS, functionDBSNP, aminoAcids, proteinPosition, distanceToSplice, polyPhen, granthamScore, scorePhastCons, scoreCADD, \
\nF_rarest_allele, MAFinESP, FILTERVCF, \
\nGenotype_QualityVCF, Variant_Confidence_Quality_by_DepthVCF, \
\nDepth_CountedVCF depth_proband, Reference_readsVCF Ref_reads_proband, Alt_readsVCF Alt_read_proband, \
\nRef_read_percentVCF Ref_read_pct_proband, A_readVCF A_read_proband, G_readVCF G_read_proband, \
\nT_readVCF T_read_proband, C_readVCF C_read_proband, A_percentVCF A_percent_proband, G_percentVCF G_percent_proband, \
\nT_percentVCF T_percent_proband, C_percentVCF C_percent_proband from " + self.proband + self.group + "; \n"  
        
        query += "ALTER TABLE " + self.proband + self.group + "_genotemp1 ADD INDEX (uniLoc), ADD INDEX (geneList), ADD INDEX(varMD5);\n"
        
        query += "CREATE TABLE " + self.proband + self.group + "_genotypes as select " + self.proband + self.group + "_genotemp1.*, \
\n" + self.parent1 + "VCF.GenotypeVCF p1_geno, " + self.parent1 + "VCF.Depth_CountedVCF Depth_p1, \
\n" + self.parent1 + "VCF.Reference_readsVCF Ref_read_p1, " + self.parent1 + "VCF.Alt_readsVCF Alt_read_p1, \
\n" + self.parent1 + "VCF.Ref_read_percentVCF Ref_read_pct_p1, " + self.parent1 + "VCF.A_readVCF A_read_p1, \
\n" + self.parent1 + "VCF.G_readVCF G_read_p1, " + self.parent1 + "VCF.T_readVCF T_read_p1, " + self.parent1 + "VCF.C_readVCF C_read_p1,\
\n" + self.parent1 + "VCF.A_percentVCF A_percent_p1, " + self.parent1 + "VCF.G_percentVCF G_percent_p1, " + self.parent1 + "VCF.T_percentVCF T_percent_p1, \
\n" + self.parent1 + "VCF.C_percentVCF C_percent_p1 from " + self.proband + self.group + "_genotemp1 left join " + self.parent1 + "VCF on \
\n" + self.proband + self.group + "_genotemp1.varMD5 = " + self.parent1 + "VCF.varMD5VCF;\n"
        
        query += "DROP TABLE " + self.proband + self.group + "_genotemp1;\n"
        
        query += "ALTER TABLE " + self.proband + group  + "_genotypes ADD INDEX (uniLoc), ADD INDEX (geneList), ADD INDEX(varMD5);\n"
        
        self.genotable = self.proband + self.group + "_genotypes"
        
        return query
    
    def getVariants(self, stringent = False):
        query = ""
        query += "SELECT * FROM " + self.proband + self.group + "_genotypes where \
\n((substring_index(proband_geno, '/', 1) = substring_index(proband_geno, '/', -1)) and \
\n(substring_index(p1_geno, '/', 1) != substring_index(p1_geno, '/', -1))) #loss of heterozygosity \
\nor ((substring_index(proband_geno, '/', 1) != substring_index(proband_geno, '/', -1)) and \
\n(substring_index(p1_geno, '/', 1) = substring_index(p1_geno, '/', -1))) #gain of heterozygosity \
\nor (substring_index(proband_geno, '/', 1) not in \
\n(substring_index(p1_geno, '/', 1), substring_index(p1_geno, '/', -1))) #one allele not found in the parent line\
\nor (substring_index(proband_geno, '/', -1) not in \
\n(substring_index(p1_geno, '/', 1), substring_index(p1_geno, '/', -1))) #the other allele not found in the parent line"
        if not stringent:
            query += ";\n"
        if stringent:
            query += "\nand proband_geno != './.' and p1_geno != './.'; #requires that genotypes were called in both parent and proband;\n;"
            
        return query
    
    def countVariants(self, stringent = False, unique = False):
        query = ""
        if unique:
            query += "SELECT count(distinct varMD5) FROM "
        if not unique:
            query += "SELECT count(*) FROM "
        query += self.proband + self.group + "_genotypes where \
\n((substring_index(proband_geno, '/', 1) = substring_index(proband_geno, '/', -1)) and \
\n(substring_index(p1_geno, '/', 1) != substring_index(p1_geno, '/', -1))) #loss of heterozygosity \
\nor ((substring_index(proband_geno, '/', 1) != substring_index(proband_geno, '/', -1)) and \
\n(substring_index(p1_geno, '/', 1) = substring_index(p1_geno, '/', -1))) #gain of heterozygosity \
\nor (substring_index(proband_geno, '/', 1) not in \
\n(substring_index(p1_geno, '/', 1), substring_index(p1_geno, '/', -1))) #one allele not found in the parent line\
\nor (substring_index(proband_geno, '/', -1) not in \
\n(substring_index(p1_geno, '/', 1), substring_index(p1_geno, '/', -1))) #the other allele not found in the parent line"
        if not stringent:
            query += ";\n"
        if stringent:
            query += "\nand proband_geno != './.' and p1_geno != './.'; #requires that genotypes were called in both parent and proband;\n;"
            
        return query
        


#code below is just for testing code above
group = "_rare_damage"
analysis = findChanges("as022plus","as022minus","_rare_damage")
print(analysis.build())
print(analysis.run(True))
print("someting")