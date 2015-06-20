#!/usr/bin/env python

class KnownGenes(object):
    
    def __init__(self, table, knownGeneList):
        self.table = table
        try:
            assert hasattr(knownGeneList, '__iter__') and not type(knownGeneList) is str
        except:
            raise TypeError("List of known genes must be passed as a list or a tupple")
        self.knownGeneList = list(knownGeneList)
        self.knownGeneList = [gene.strip() for gene in self.knownGeneList]
        self.knownGeneList = ["'" + gene + "'," for gene in self.knownGeneList]
        
        
        
    def listFromFile(file):
        if not hasattr(file, 'read'):
            raise TypeError("Object passed must be a file")
        linecollector = []
        line = file.readline()
        while line:
            if not line[0] == '#':
                linecollector.append(line.strip())
            line = file.readline()
        return linecollector
        
        
    def getVariants(self):
        knownGeneString = " ".join(self.knownGeneList)
        knownGeneString = knownGeneString.strip()
        knownGeneString = knownGeneString.strip(",")
        query = ""
        query += "SELECT * FROM " + self.table + " where geneList in (" + knownGeneString + ") ORDER by geneList;"
        return query
    
    
patient = "patient"
file = open('skeletalgenelist.txt', 'r')
genelist = KnownGenes.listFromFile(file)
query = KnownGenes(patient, genelist)
print(query.getVariants())
print(str(len(query.knownGeneList)))
print("something")