#!/usr/bin/env python

class Homozygous(object):
    
    def __init__(self, table, readCutoff = 0.1, suffix = "_homozygosity"):
        self.inputTable = table
        self.tablename = table + suffix
        if readCutoff:
            try:
                readCutoff = float(readCutoff)
            except TypeError:
                raise TypeError("Read cutoff value must be a number.")
            if readCutoff >= 1.0:
                readCutoff = readCutoff/100
            try:
                assert readCutoff <= 1.0
            except:
                raise ValueError("Percent for cutoff must be below 100")
            if readCutoff > 0.5:
                readCutoff = 1.0 - readCutoff
        self.readCutoff = readCutoff
        
    def buildTables(self):
        query = ""
        query += "CREATE TABLE " + self.tablename + " as SELECT * from " + self.inputTable + " \
\nWHERE sampleGenotype NOT IN ('S' , 'R', 'Y', 'M', 'K', 'W') \
\nAND (INSTR(sampleGenotype, '/') = 0 OR substring_index(sampleGenotype, '/', 1) = substring_index(sampleGenotype, '/', -1))"
        if self.readCutoff:
            query += "\nOR (Ref_read_percentVCF < " + str(round(self.readCutoff, 4)) + ");\n"
        else:
            query += ";\n"
        query += "ALTER TABLE " + self.tablename + " ADD INDEX (varMD5), ADD INDEX (uniLoc), ADD INDEX (geneList);"
        
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
    
#begin testcode below
patient = "patient"
query = Homozygous(patient, 95)
print(query.buildTables())
print(query.countVariants())
print(query.tablename)
print("something")
        