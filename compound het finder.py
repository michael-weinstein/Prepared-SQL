#!/usr/bin/env python

class CompoundHet(object):
    
    def __init__(self, table, suffix = "_compound_het"):
        self.inputTable = table
        self.tablename = table + suffix
        
    def buildTables(self):
        hitcounts = self.inputTable + "_hit_count"
        query = ""
        query += "CREATE TABLE " + hitcounts + " as \
\nSELECT distict varMD5, uniLoc, geneList, count(*) hits \
\nFROM " + self.inputTable + " group by geneList; \
\nCREATE TABLE " + self.tablename + " as SELECT * from " + self.inputTable + " \
\nWHERE geneList in (SELECT distinct geneList from " + hitcounts + " WHERE \
hits > 1) \
\nORDER by geneList asc;"
        
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
query = CompoundHet(patient)
print(query.buildTables())
print(query.countVariants())
print(query.tablename)
print("something")