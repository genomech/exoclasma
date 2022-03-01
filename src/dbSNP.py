from .SharedFunctions import * 
from .TabixFunctions import *


# dbSNP v151

dbSNPFlags = ["RV", "PM", "TPA", "PMC", "S3D", "SLO", "NSF", "NSM", "NSN", "REF", "SYN", "U3",
			  "U5", "ASS", "DSS", "INT", "R3", "R5", "OTH", "CFL", "ASP", "MUT", "VLD", "G5A",
			  "G5", "HD", "GNO", "KGPhase1", "KGPhase3", "CDA", "LSD", "MTP", "OM", "NOC", "WTD", "NOV"]

def dbSNPParseLine(Row):
	Result = ParseVcfRow(Row)
	NewInfo = {}
	Keys = list(Result['INFO'].keys())
	Alleles = [Result['REF']] + Result['ALT']
	Flags = []
	for Key in Keys:
		if Key == 'RS':           NewInfo['RS']         = int(Result['INFO']['RS'])
		if Key == 'RSPOS':        NewInfo['RSPOS']      = int(Result['INFO']['RSPOS'])
		if Key == 'VP':           NewInfo['VP']         = str(Result['INFO']['VP'])
		if Key == 'GENEINFO':     NewInfo['GENEINFO']   = [{ 'Gene': i.split(':')[0], 'ID': i.split(':')[1] } for i in Result['INFO']['GENEINFO'].split('|')]
		if Key == 'dbSNPBuildID': NewInfo['dbSNPBuild'] = int(Result['INFO']['dbSNPBuildID'])
		if Key == 'SAO':          NewInfo['SAO']        = int(Result['INFO']['SAO'])
		if Key == 'SSR':          NewInfo['SSR']        = int(Result['INFO']['SSR'])
		if Key == 'WGT':          NewInfo['WGT']        = int(Result['INFO']['WGT'])
		if Key == 'VC':           NewInfo['VC']         = str(Result['INFO']['VC'])
		if Key == 'CAF':          NewInfo['CAF']        = { Alleles[index]: None if (item == '.') else float(item) for index, item in enumerate(Result['INFO']['CAF'].split(',')) }
		if Key == 'TOPMED':       NewInfo['TOPMED']     = { Alleles[index]: None if (item == '.') else float(item) for index, item in enumerate(Result['INFO']['TOPMED'].split(',')) }
		
		if Key in dbSNPFlags: Flags += [Key]
		if Key == 'COMMON':
			if Result['INFO'][Key] == 1: Flags += [Key]
	NewInfo['Flags'] = Flags
	Result['INFO'] = NewInfo
	return Result

def dbSNPVariantMatch(Item, Row):
	return (Item[0] == Row[0]) and (int(Item[1]) == int(Row[1])) and (Item[2] == Row[3]) and (Item[3] in Row[4].split(','))

def dbSNPQueryFunction(Item, TabixObj):
	TabixOutput = TabixObj.query(Item[0], Item[1] - 1, Item[1])
	Result = [ dbSNPParseLine(Row = Row) for Row in TabixOutput if dbSNPVariantMatch(Item = Item, Row = Row) ]
	BriefDetails = { "dbSNP": None if not Result else [item['ID'] for item in Result], "dbSNP_Details": None if not Result else Result }
	return Result
