#!/usr/bin/env python3

class RareCodingVariants(object):
    
    def __init__(self, table, suffix = "_rare_coding"):
        self.inputTable = table
        self.tablename = table + suffix
        
    def buildTables(self):
        temptable = self.tablename + "_temp"
        query = ""
        query += "CREATE TABLE " + temptable + " as SELECT * from " + self.inputTable + " where \
\n(AfricanHapMapFreq < 1.0) and (AsianHapMapFreq < 1.0) and (EuropeanHapMapFreq < 1.0) and \
\n(dbSNP_5_percent_in_all != 'TRUE') and (dbSNP_5_percent_in_any != 'TRUE') and \
\n(F_rarest_allele < 0.01) and (MAFinESP < 0.01) and (functionGVS NOT IN \
\n('intron' , 'synonymous', '5-prime-UTR', '3-prime-UTR', 'non-coding-exon', \
\n'upstream-gene', 'downstream-gene', 'intergenic')) \
\nGROUP by varMD5; \
\nDELETE from " + self.inputTable + "VCF where instr(uniLocVCF, 'GL') != 0; \
\nCREATE TABLE " + self.tablename + " as select * from " + temptable + " left join " + self.inputTable + "VCF on varMD5 = varMD5VCF; \
\nDROP TABLE " + temptable + "; \
\nALTER TABLE " + self.tablename + " ADD INDEX (varMD5), ADD INDEX (uniLoc), ADD INDEX (geneList);"
        
        return query
    
    
    def getVariants(self):
        query = "SELECT * FROM " + self.tablename + ";"
        return query
    
    
    def countVariants(self, unique = False):
        query = ""
        if unique:
            query += "SELECT COUNT(distinct varMD5) "
        if not unique:
            query += "SELECT COUNT(*) "
        query += "from " + self.tablename + ";"
        
        return query
        


#testcode below
patient = "patient"
query = RareCodingVariants(patient)
print(query.buildTables())
print(query.countVariants())
print(query.tablename)
print("something")