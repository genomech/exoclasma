from .SharedFunctions import * 
from .TabixFunctions import *


# ClinVar 2022/02/27

def ClinVarParseLine(Row):
	Result = ParseVcfRow(Row)
	NewInfo = {}
	Keys = list(Result['INFO'].keys())
	for Key in Keys:
		if Key == 'AF_ESP': NewInfo['AF_ESP'] = float(Result['INFO']['AF_ESP'])
		if Key == 'AF_EXAC': NewInfo['AF_EXAC'] = float(Result['INFO']['AF_EXAC'])
		if Key == 'AF_TGP': NewInfo['AF_TGP'] = float(Result['INFO']['AF_TGP'])
		if Key == 'ALLELEID': NewInfo['ALLELEID'] = int(Result['INFO']['ALLELEID'])
		if Key == 'CLNDN': NewInfo['CLNDN'] = [item.replace('_', ' ') for item in Result['INFO']['CLNDN'].split('|')]
		if Key == 'CLNDNINCL': NewInfo['CLNDNINCL'] = [item.replace('_', ' ') for item in Result['INFO']['CLNDNINCL'].split('|')]
		if Key == 'CLNDISDB': NewInfo['CLNDISDB'] = None if set(Result['INFO']['CLNDISDB'].split('|')) == {'.'} else [ {'DB': item.split(':', maxsplit=1)[0], 'ID': item.split(':', maxsplit=1)[1]} for item in Result['INFO']['CLNDISDB'].split(',') ]
		if Key == 'CLNDISDBINCL': NewInfo['CLNDISDBINCL'] = None if set(Result['INFO']['CLNDISDBINCL'].split('|')) == {'.'} else [ {'DB': item.split(':', maxsplit=1)[0], 'ID': item.split(':', maxsplit=1)[1]} for item in Result['INFO']['CLNDISDBINCL'].split(',') ]
		if Key == 'CLNHGVS': NewInfo['CLNHGVS'] = str(Result['INFO']['CLNHGVS'])
		if Key == 'CLNREVSTAT': NewInfo['CLNREVSTAT'] = Result['INFO']['CLNREVSTAT'].replace('_', ' ')
		if Key == 'CLNSIG': NewInfo['CLNSIG'] = Result['INFO']['CLNSIG'].replace('_', ' ')
		if Key == 'CLNSIGCONF': NewInfo['CLNSIGCONF'] = [ { "ClnSig": item.split('(')[0].replace('_', ' '), "Count": int(item.split('(')[1].split(')')[0]) } for item in Result['INFO']['CLNSIGCONF'].split(',') ]
		if Key == 'CLNSIGINCL': NewInfo['CLNSIGINCL'] = [ {'VariationID': item.split(':', maxsplit=1)[0], 'ClnSig': item.split(':', maxsplit=1)[1] } for item in Result['INFO']['CLNSIGINCL'].split('|') ]
		if Key == 'CLNVC': NewInfo['CLNVC'] = Result['INFO']['CLNVC'].replace('_', ' ')
		if Key == 'CLNVCSO': NewInfo['CLNVCSO'] = Result['INFO']['CLNVCSO']
		if Key == 'CLNVI': NewInfo['CLNVI'] = [{'DB': item.split(':', maxsplit=1)[0].replace('_', ' '), 'ID': item.split(':', maxsplit=1)[1]} for item in Result['INFO']['CLNVI'].split('|')]
		if Key == 'DBVARID': NewInfo['DBVARID'] = Result['INFO']['DBVARID']
		if Key == 'GENEINFO': NewInfo['GENEINFO'] = [{ 'Gene': i.split(':')[0], 'ID': i.split(':')[1] } for i in Result['INFO']['GENEINFO'].split('|')]
		if Key == 'MC': NewInfo['MC'] = [{"ID": item.split('|')[0], "Consequence": item.split('|')[1].replace('_', ' ')} for item in Result['INFO']['MC'].split(',')]
		if Key == 'ORIGIN': NewInfo['ORIGIN'] = int(Result['INFO']['ORIGIN'])
		if Key == 'RS': NewInfo['RS'] = [int(item) for item in Result['INFO']['RS'].split('|')]
		if Key == 'SSR': NewInfo['SSR'] = int(Result['INFO']['SSR'])
	Result['INFO'] = NewInfo
	return Result

def ClinVarQueryFunction(Item, TabixObj):
	TabixOutput = TabixObj.query(Item[0], Item[1] - 1, Item[1])
	Result = [ ClinVarParseLine(Row = Row) for Row in TabixOutput if VcfVariantMatch(Item = Item, Row = Row) ]
	return None if not Result else Result

