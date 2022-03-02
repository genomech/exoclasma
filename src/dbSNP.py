from .SharedFunctions import * 
from .TabixFunctions import *


# dbSNP v151

dbSNPFlags = ["RV", "PM", "TPA", "PMC", "S3D", "SLO", "NSF", "NSM", "NSN", "REF", "SYN", "U3",
			 "U5", "ASS", "DSS", "INT", "R3", "R5", "OTH", "CFL", "ASP", "MUT", "VLD", "G5A",
			 "G5", "HD", "GNO", "KGPhase1", "KGPhase3", "CDA", "LSD", "MTP", "OM", "NOC", "WTD", "NOV"]
dbSNPSimpleInt = ['RS', 'RSPOS', 'dbSNPBuildID', 'SAO', 'SSR', 'WGT']

def dbSNPParseLine(Row):
	Result = ParseVcfRow(Row)
	NewInfo = {}
	Keys = list(Result['INFO'].keys())
	Alleles = [Result['REF']] + Result['ALT']
	Flags = []
	for Key in Keys:
		if Key == 'VP': NewInfo['VP'] = str(Result['INFO']['VP'])
		if Key == 'GENEINFO': NewInfo['GENEINFO'] = [{ 'Gene': i.split(':')[0], 'ID': i.split(':')[1] } for i in Result['INFO']['GENEINFO'].split('|')]
		if Key == 'VC': NewInfo['VC'] = str(Result['INFO']['VC'])
		if Key == 'CAF': NewInfo['CAF'] = { Alleles[index]: None if (item == '.') else float(item) for index, item in enumerate(Result['INFO']['CAF'].split(',')) }
		if Key == 'TOPMED': NewInfo['TOPMED'] = { Alleles[index]: None if (item == '.') else float(item) for index, item in enumerate(Result['INFO']['TOPMED'].split(',')) }
		if Key in dbSNPSimpleInt: NewInfo[Key] = int(Result['INFO'][Key])
		if Key in dbSNPFlags: Flags += [Key]
		if Key == 'COMMON':
			if Result['INFO'][Key] == 1: Flags += [Key]
	NewInfo['Flags'] = Flags
	Result['INFO'] = NewInfo
	return Result

def dbSNPQueryFunction(Item, TabixObj):
	TabixOutput = TabixObj.query(Item[0], Item[1] - 1, Item[1])
	Result = [ dbSNPParseLine(Row = Row) for Row in TabixOutput if VcfVariantMatch(Item = Item, Row = Row) ]
	return None if not Result else Result
