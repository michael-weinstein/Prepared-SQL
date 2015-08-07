#!/usr/bin/env python3
#query building functions and objects
class MapReduce(object):
    
    def __init__(self, data):
        self.data = [str(datum) for datum in data]
        self.generate()
        
    def generate(self):
        self.hashTable = {}
        for datum in self.data:
            try:
                self.hashTable[datum] += 1
            except KeyError:
                self.hashTable[datum] = 1
        keyList = list(self.hashTable.keys())
        keyList.sort()
        self.map = keyList
        self.reduced = []
        for key in self.map:
            self.reduced.append([key,self.hashTable[key]])
        return True

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

    
class KnownGenes(object):
    
    def __init__(self, table, knownGeneList, columnName = 'geneList'):
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
    
    def hashedListFromFile(file):
        if not hasattr(file, 'read'):
            raise TypeError("Object passed must be a file")
        listcollector = {}
        line = file.readline()
        line = line.strip()
        if line[0] != "#":
            raise RuntimeError("First line must be a group name and start with a # symbol")
        while line:
            if line[0] == "#":
                currentlist = line.replace("#","")
                currentlist = currentlist.replace("\t"," ")
                listcollector[currentlist] = []
            else:
                listcollector[currentlist].append(line.strip())
            line = file.readline()
            line = line.strip()
        return listcollector
       
    def getVariants(self, columnName = "geneList"):
        knownGeneString = " ".join(self.knownGeneList)
        knownGeneString = knownGeneString.strip()
        knownGeneString = knownGeneString.strip(",")
        query = ""
        query += "SELECT * FROM " + self.table + " where " + columnName + " in (" + knownGeneString + ") ORDER by geneList;"
        return query
    

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
    
    def countGenes(self):
        query = ""
        query += "SELECT COUNT(distinct geneList)from " + self.tablename + ";"
        return query
    
    
class DeNovoChanges(object):
    
    def __init__(self, proband, parent1 = "", parent2 = "", group = ""):
        if not (proband and parent1 and parent2) or not type(proband) is str:
            try:
                assert hasattr(proband, '__iter__') and len(proband) == 3
                self.group = parent1
                self.parent2 = proband[2]
                self.parent1 = proband[1]
                self.proband = proband[0]
            except:
                raise RuntimeError("Finding changed variants requires passing proband and original as either arguments or a tupple of two identifiers")
        else:
            self.proband = proband
            self.group = group
            self.parent1 = parent1
            self.parent2 = parent2            
        
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
        query += "CREATE TABLE " + self.proband + self.group + "_genotemp2 as select " + self.proband + self.group + "_genotemp1.*, \
\n" + self.parent1 + "VCF.GenotypeVCF p1_geno, " + self.parent1 + "VCF.Depth_CountedVCF Depth_p1, \
\n" + self.parent1 + "VCF.Reference_readsVCF Ref_read_p1, " + self.parent1 + "VCF.Alt_readsVCF Alt_read_p1, \
\n" + self.parent1 + "VCF.Ref_read_percentVCF Ref_read_pct_p1, " + self.parent1 + "VCF.A_readVCF A_read_p1, \
\n" + self.parent1 + "VCF.G_readVCF G_read_p1, " + self.parent1 + "VCF.T_readVCF T_read_p1, " + self.parent1 + "VCF.C_readVCF C_read_p1,\
\n" + self.parent1 + "VCF.A_percentVCF A_percent_p1, " + self.parent1 + "VCF.G_percentVCF G_percent_p1, " + self.parent1 + "VCF.T_percentVCF T_percent_p1, \
\n" + self.parent1 + "VCF.C_percentVCF C_percent_p1 from " + self.proband + self.group + "_genotemp1 left join " + self.parent1 + "VCF on \
\n" + self.proband + self.group + "_genotemp1.varMD5 = " + self.parent1 + "VCF.varMD5VCF;\n"
        query += "DROP TABLE " + self.proband + self.group + "_genotemp1;\n"
        query += "ALTER TABLE " + self.proband + group  + "_genotemp2 ADD INDEX (uniLoc), ADD INDEX (geneList), ADD INDEX(varMD5);\n"
        query += "CREATE TABLE " + self.proband + self.group + "_genotypes as select " + self.proband + self.group + "_genotemp2.*, \
\n" + self.parent2 + "VCF.GenotypeVCF p2_geno, " + self.parent2 + "VCF.Depth_CountedVCF Depth_p2, \
\n" + self.parent2 + "VCF.Reference_readsVCF Ref_read_p2, " + self.parent2 + "VCF.Alt_readsVCF Alt_read_p2, \
\n" + self.parent2 + "VCF.Ref_read_percentVCF Ref_read_pct_p2, " + self.parent2 + "VCF.A_readVCF A_read_p2, \
\n" + self.parent2 + "VCF.G_readVCF G_read_p2, " + self.parent2 + "VCF.T_readVCF T_read_p2, " + self.parent2 + "VCF.C_readVCF C_read_p2,\
\n" + self.parent2 + "VCF.A_percentVCF A_percent_p2, " + self.parent2 + "VCF.G_percentVCF G_percent_p2, " + self.parent2 + "VCF.T_percentVCF T_percent_p2, \
\n" + self.parent2 + "VCF.C_percentVCF C_percent_p2 from " + self.proband + self.group + "_genotemp2 left join " + self.parent2 + "VCF on \
\n" + self.proband + self.group + "_genotemp2.varMD5 = " + self.parent2 + "VCF.varMD5VCF;\n"
        query += "DROP TABLE " + self.proband + self.group + "_genotemp2;\n"
        query += "ALTER TABLE " + self.proband + group  + "_genotypes ADD INDEX (uniLoc), ADD INDEX (geneList), ADD INDEX(varMD5);\n"
        self.genotable = self.proband + self.group + "_genotypes"
        return query
    
    def getVariants(self, stringent = False):
        query = ""
        query += "SELECT * FROM " + self.proband + self.group + "_genotypes WHERE \
\n(substring_index(proband_geno, '/', 1) != substring_index(proband_geno, '/', -1)) and \
\n((substring_index(proband_geno, '/', 1) != '0' and substring_index(proband_geno, '/', 1) not in \
\n (substring_index(mom_geno, '/', 1), substring_index(mom_geno, '/', -1), \
\n  substring_index(dad_geno, '/', 1), substring_index(dad_geno, '/', -1))) \
\nor \
\n(substring_index(proband_geno, '/', -1) != '0' and substring_index(proband_geno, '/', -1) not in \
\n (substring_index(mom_geno, '/', 1), substring_index(mom_geno, '/', -1), \
\n  substring_index(dad_geno, '/', 1), substring_index(dad_geno, '/', -1)) \
\n))"
        if not stringent:
            query += ";\n"
        if stringent:
            query += "\nand proband_geno != './.' and p1_geno != './.' and p2_geno != './.';\n"
        return query   
    
    def countVariants(self, stringent = False, unique = False):
        query = ""
        if unique:
            query += "SELECT COUNT(distinct varMD5) FROM "
        if not unique:
            query += "SELECT COUNT(*) FROM "
        query += self.proband + self.group + "_genotypes WHERE \
\n(substring_index(proband_geno, '/', 1) != substring_index(proband_geno, '/', -1)) and \
\n((substring_index(proband_geno, '/', 1) != '0' and substring_index(proband_geno, '/', 1) not in \
\n (substring_index(mom_geno, '/', 1), substring_index(mom_geno, '/', -1), \
\n  substring_index(dad_geno, '/', 1), substring_index(dad_geno, '/', -1))) \
\nor \
\n(substring_index(proband_geno, '/', -1) != '0' and substring_index(proband_geno, '/', -1) not in \
\n (substring_index(mom_geno, '/', 1), substring_index(mom_geno, '/', -1), \
\n  substring_index(dad_geno, '/', 1), substring_index(dad_geno, '/', -1)) \
\n))"
        if not stringent:
            query += ";\n"
        if stringent:
            query += "\nand proband_geno != './.' and p1_geno != './.' and p2_geno != './.';\n"
        return query 
    
    
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
    
    
class CleanUp(object):
    
    def rowstring(data, delimiter = "\t", prepend = "", append = ""):
        rowcollector = ""
        for row in range(0,len(data)):
            if delimiter == "\t":
                if prepend:
                    columncollector = prepend + delimiter
                else:
                    columncollector = ""
            else:
                if prepend:
                    columncollector = '"' + prepend + '"' + delimiter + '"'
                else:
                    columncollector = '"'
            for column in range(0,len(data[row])):
                columncollector += str(data[row][column])
                if not (len(data) == 1 and len(data[row]) == 1):
                    if column == len(data[row]) -1:
                        if delimiter == "\t":
                            if append:
                                columncollector += delimiter + append + delimiter + "\n"
                            else:
                                columncollector += "\n"
                        else:
                            if append:
                                columncollector += '"' + delimiter + '"' + append + '"' + "\n"
                            else:
                                columncollector += '"' + "\n"
                    else:
                        if delimiter == "\t":
                            columncollector += delimiter
                        else:
                            columncollector += '"' + delimiter + '"'
            rowcollector += columncollector
        return rowcollector
    
    def columnstring(data, delimiter = "\t"):
        if delimiter == "\t":
            rowcollector = ""
        else:
            rowcollector = '"'
        for row in range(0,len(data)):
            rowcollector += "".join(data[row])
            if row == len(data) -1:
                if delimiter == "\t":
                    rowcollector += "\n"
                else:
                    rowcollector += '"' + "\n"
            else:
                if delimiter == "\t":
                    rowcollector += delimiter
                else:
                    rowcollector += '"' + delimiter + '"'
        return rowcollector
    
    def returnToArray(data, singleColumn = False):
        rowcollector = []
        for row in data:
            columncollector = []
            for column in row:
                columncollector.append(column)
            rowcollector.append(columncollector)
        if singleColumn:
            rowcollector = [item[0] for item in rowcollector]
        return rowcollector
        
        
class Sample(object):
    
    def __init__(self, samples):
        try:
            self.samples = int(samples)
        except:
            raise ValueError("Number of samples to collect must be an integer (or string of an integer, if you really must)")
        self.collected = 0
        self.collection = []
        
    def collect(self, data):
        if self.collected < self.samples:
            self.collection.append(data)
            self.collected += 1
            
    def display(self):
        for item in self.collection:
            print(item)
            
            
class JoinWithVCFTable(object):
    
    def __init__(self, table, suffix = "_all"):
        self.table = table
        self.suffix = suffix
        
    def buildTables(self):
        query = ""
        query += "CREATE TABLE " + self.table + self.suffix + " as \
\nSELECT * from " + self.table + " left join " + self.table + "VCF on varMD5 = varMD5VCF;"
        self.joinedTable = self.table + self.suffix
        return query
    
def getColumns(schema, table):
    query = "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='" + schema + "' AND `TABLE_NAME`='" + table + "';"
    return query
    
def queryBuilderHelp():
    print("RareCodingVariants(tablename, suffix) Finds variants annotated as less than 1% allele frequency in all tested populations that are annotated as variant types with the potential to affect protein sequence. Suffix for new table defaults to _rare_coding if not set\n\tbuildTables()\n\tgetVariants()\n\tcountVariants()")
    print("KnownGenes(tablename, knownGeneList) Finds variants from the given table in the list of known genes.  knownGeneList should be a list or tupple.\n\tlistFromFile(fileObject) (builds a list of genes from a text file, ignores lines starting with a # symbol)\n\thashedListFromFile(fileObject) (builds multiple lists of files under one dictionary using lines starting with a # symbol as the key)\n\tgetVariants()")
    print("Homozygous(tablename, readCutoff, suffix) Creates a table of homozygous variants from the tablename entered.\nRead cutoff allows for including variants not called as homozygous by the variant caller.\nNew table suffix defaults to _homozygosity.\n\tbuildTables()\n\tgetVariants()\n\tcountVariants(unique) If the same variant is listed multiple times due to listing each transcript, setting this to True will only count it once (default value is False)")
    print("CompoundHet(tablename, suffix) Creates a table of genes with multiple variants.  Default suffix is _compound_het.  VARIANTS MAY OR MAY NOT BE IN PHASE\n\tbuildTables()\n\tgetVariants()\n\tcountVariants(unique)\n\tcountGenes() Returns a count of the genes with multiple variants")
    print("DeNovoChanges(proband, parent1, parent2, group) Finds variants that have a variant genotype in the proband not found in one of the parents.  Group defaults to an empty string, but can be used to add a suffix to the proband's tablename.\n\tbuildTables()\n\tgetVariants(stringent) Stringent is a boolean value (True or False) with a default value of False.  If true, will leave out any loci that were not called in one or both parents.\n\tcountVariants(stringent, unique)")
    print("SomaticChanges(proband, parent1, group) Similar to above, but looks between proband and one parent for any genotypes in the proband not found in the parent, gains of heterozygosity, or losses of heterozygosity.\n\tbuildTables()\n\tgetVariants(stringent)\n\tcountVariants(stringent, unique)")
    print("CleanUp() A holder for functions that turn the SQL output into clean strings or lists\n\trowstring(data, delimiter) Returns a delimited string (default is tab) of the submitted data.\n\tcolumnstring(data, delimiter) Turns a 2-dimensional output into a 1-line dimensional string (useful for column headers)\n\trowarray(data) Returns a 2-dimensional array of data\n\tcolumnarray(data)")
    print("Sample(numberOfSamples) Sets up an object to store the number of samples given.\n\tcollect(data) Stores the data as a sample.\n\tdisplay() Prints the collected samples.")
    print("JoinWithVCFTable(tablename, suffix) Joins the annotation (SeattleSeq) data with the matching VCF data for that variant in the tablename submitted.  Suffix default is _all.\n\tbuildTables()")
    print("getColumns(schema, tablename) Gets a list of column names from the schema and tablename.  Data returned as a column of tupples each containing one column name. (Use columnstring or column array to clean).")   


class Recurrents(object):
    
    def __init__(self, caseDB, caseTableList, controlDB=False, controlTableList=[]):
        self.cases = caseTableList
        self.controls = controlTableList
        self.caseDB = caseDB
        self.controlDB = controlDB
        self.caseCounts = self.getCaseCounts(self.caseDB, self.cases)
        if self.controlDB and self.controls:
            self.controlCounts = self.getControlCounts(self.controlDB, self.controls)
        else:
            self.controlCounts = False
            
    def getCaseCounts(self, dbConnection, tables):
        md5collector = []
        for table in tables:
            query = "Select varMD5 from " + table
            dbConnection.execute(query)
            md5list = dbConnection.fetchall()
            md5list = CleanUp.returnToArray(md5list, "SingleColumn")
            md5collector += md5list
        counts = MapReduce(md5collector)
        self.caseRecurrenceList = counts.reduced
        self.caseRecurrenceHash = counts.hashTable
        self.caseVariantList = counts.map
        
    def getControlCounts(self, dbConnection, tables):
        md5collector = []
        for table in tables:
            query = "Select varMD5 from " + table
            dbConnection.execute(query)
            md5list = dbConnection.fetchall()
            md5list = CleanUp.returnToArray(md5list, "SingleColumn")
            md5collector += md5list
        counts = MapReduce(md5collector)
        self.controlRecurrenceList = counts.reduced
        self.controlRecurrenceHash = counts.hashTable
        self.controlVariantList = counts.map
        
    def hashByCount(self, control = False):
        output = {}
        if control:
            hashTable = self.controlRecurrenceHash
            maxRecurrences = len(self.controls)
        else:
            hashTable = self.caseRecurrenceHash
            maxRecurrences = len(self.cases)
        for count in range(1,maxRecurrences + 1):
            output[count] = []
        for key in list(hashTable.keys()):
            output[hashTable[key]].append(key)
        return output
    
class multiTableFind(object):
    def __init__(self, dbconnection, schema, columnName, matchCriteria, tableList, returnColumns = "*"):
        self.columnName = columnName
        self.matchCriteria = matchCriteria
        self.matchCriteriaString = "'" + "', '".join(matchCriteria) + "'"
        self.tableList = tableList
        self.returnColumns = returnColumns
        self.dbconnection = dbconnection
        self.schema = schema
        self.returnColumnsString = "*"
        if self.returnColumns != "*":
            self.returnColumnsString = ", ".join(returnColumns)
            self.headerList = returnColumns
        if not self.validateColumns():
            if self.returnColumns != "*":
                raise RuntimeError("The following requested columns are missing (shown as table.column):\n" + "\n".join(self.missingColumns))
            else:
                tableDebug = ""
                for table in self.tableList:
                    tableDebug += "\t".join([table,self.tableDebugHash[table][0],self.tableDebugHash[table][1]]) + "\n"
                raise RuntimeError("All tables must have the same headers if all columns are requested in the return.\
                \nDebugging information as follows:\n" + tableDebug)
        else:
            self.run()
                
    def validateColumns(self):
        import hashlib
        passedCheck = True
        if self.returnColumns == "*":
            self.tableDebugHash = {}
            firstHeader = []
            for table in self.tableList:
                self.dbconnection.execute(getColumns(self.schema, table))
                result = self.dbconnection.fetchall()
                result = CleanUp.returnToArray(result, singleColumn = True)
                result[88] = "original_sample_info"  #if annotation columns change, this will have to be changed with it to a new index
                headerString = "".join(result)
                self.tableDebugHash[table] = [str(len(result)), hashlib.md5(headerString.encode('utf-8')).hexdigest()]
                #print(result)
                if not firstHeader:
                    firstHeader = result
                    self.headerList = result
                else:
                    for i in range(0,len(firstHeader)):
                        if result[i] != firstHeader[i]:
                            print(str(i) + "\n" + result[i] + "\n" + firstHeader[i])
                            passedCheck = False
            #        if result != firstHeader:
            #            passedCheck = False
        else:
            self.missingColumns = []
            for table in self.tableList:
                self.dbconnection.execute(getColumns(self.schema, table))
                result = self.dbconnection.fetchall()
                result = CleanUp.returnToArray(result, singleColumn = True)
                for column in self.returnColumns:
                    if not column in result:
                        self.missingColumns.append(".".join([table,column]))
                        passedCheck = False
        return passedCheck
    
    
    def run(self):
        self.result = {}
        for table in self.tableList:
            query = "SELECT " + self.returnColumnsString + " from " + table + " WHERE " + self.columnName + " in (" + self.matchCriteriaString + ");\n"
            self.dbconnection.execute(query)
            returned = self.dbconnection.fetchall()
            returned = CleanUp.returnToArray(returned)
            self.result[table] = returned
        return True
    
    def makePandas(self):
        import pandas
        thisWillProbablyFail = (ValueError, TypeError)
        pandaHeaders = ['source'] + self.headerList
        pandaPrep = {}
        for header in pandaHeaders:
            pandaPrep[header] = []
        for key in (self.result.keys()):
            for line in self.result[key]:
                position = 1 #index to 1 because we add the source value by name
                pandaPrep['source'].append(key)
                for element in line:
                    try:
                        element = int(element)
                    except thisWillProbablyFail:
                        try:
                            element = float(element)
                        except thisWillProbablyFail:
                            pass
                    pandaPrep[pandaHeaders[position]].append(element)
                    position += 1
        df =  pandas.DataFrame(pandaPrep)
        return df