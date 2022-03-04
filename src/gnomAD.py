from .SharedFunctions import * 
from .TabixFunctions import *


# gnomAD 2.1.1

gnomAD_Scheme = LoadJSON(os.path.join(CurrentDir(), '..', 'schemes', 'gnomAD_v.2.1.1_scheme.json'))

gnomAD_Functions = {
	"int": int,
	"float": float,
	"str": str,
	"gnomAD_ABHist_Parser": lambda x: {f"{gnomAD_Scheme['AB_BinEdges'][index]}-{gnomAD_Scheme['AB_BinEdges'][index + 1]}": int(item) for index, item in enumerate(x.split('|'))},
	"gnomAD_VEP_Parser": lambda x: [{ gnomAD_Scheme['VepScheme'][index]: item for index, item in enumerate(item2.split('|')) if item != ''} for item2 in x.split(',')],
	"gnomAD_AgeHist_Parser": lambda x: {f"{gnomAD_Scheme['Ages_BinEdges'][index]}-{gnomAD_Scheme['Ages_BinEdges'][index + 1]}": {"Count": int(item), "Total": gnomAD_Scheme['Ages_TotalIndividualsByBin'][index]} for index, item in enumerate(x.split('|'))},
	"gnomAD_GQnDPHist_Parser": lambda x: {f"{gnomAD_Scheme['GQnDP_BinEdges'][index]}-{gnomAD_Scheme['GQnDP_BinEdges'][index + 1]}": int(item) for index, item in enumerate(x.split('|'))}
}

def gnomAD_ParseLine(Row):
	Result = ParseVcfRow(Row)
	Alleles = [Result['REF']] + Result['ALT']
	Keys = list(Result['INFO'].keys())
	SortedKeys = { 'Flags': [], 'Tags': [] }
	for Key in Keys:
		if Key in gnomAD_Scheme["InfoFlags"]: SortedKeys['Flags'].append(Key)
		else: SortedKeys['Tags'].append(Key)
	NewInfo = {}
	for Tag in SortedKeys['Tags']:
		NewInfo[tuple(gnomAD_Scheme["InfoTags"][Tag]["Path"])] = gnomAD_Functions[gnomAD_Scheme["InfoTags"][Tag]["Func"]](Result['INFO'][Tag])
	NewInfo1 = RecursiveDict(NewInfo)
	NewInfo1['Flags'] = SortedKeys['Flags']
	Result['INFO'] = NewInfo1
	return Result

def gnomAD_QueryFunction(Item, TabixObj):
	TabixOutput = TabixObj.query(Item[0], Item[1] - 1, Item[1])
	Result = []
	for Row in TabixOutput:
		if VcfVariantMatch(Item = Item, Row = Row):
			try:
				Result.append(gnomAD_ParseLine(Row = Row))
			except Exception as Exc:
				raise ParseLineError(f"Row: {Row}; Exception Info: {Exc}")
	return None if not Result else Result
