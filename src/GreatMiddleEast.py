from .SharedFunctions import * 
from .TabixFunctions import *

# GME Variome PanTro2 [hg19]

GME_Scheme = LoadJSON(os.path.join(CurrentDir(), '..', 'schemes', 'GME_Variome_PanTro2_hg19_scheme.json'))

def GME_CalculateFrequency(string):
	Result = {}
	if '|' not in string: # Autosome
		HomRef, Het, HomAlt = [int(item) for item in string.split(':')]
		Result["AN"] = ((HomRef + Het + HomAlt) * 2)
		Result["AC"] = (Het + (HomAlt * 2))
		Result["AF"] = 0.0 if Result["AN"] == 0 else (Result["AC"] / Result["AN"])
		Result["HOM"] = HomAlt
		Result["HET"] = Het
	else: # XY
		Result = {}
		Female, Male = [item[2:-1] for item in string.split('|')]
		FemaleHomRef, FemaleHet, FemaleHomAlt = [int(item) for item in Female.split(':')]
		MaleHemRef, MaleHemAlt = [int(item) for item in Male.split(':')]
		Result["AN"] = ((FemaleHomRef + FemaleHet + FemaleHomAlt) * 2) + (MaleHemRef + MaleHemAlt)
		Result["AC"] = (FemaleHomAlt * 2) + FemaleHet + MaleHemAlt
		Result["AF"] = 0.0 if Result["AN"] == 0 else (Result["AC"] / Result["AN"])
		Result["HET"] = FemaleHet
		Result["HOM"] = FemaleHomAlt
		Result["HEM"] = MaleHemAlt
		Result["MNUM"] = MaleHemRef + MaleHemAlt
		Result["FNUM"] = FemaleHomRef + FemaleHet + FemaleHomAlt
	return Result

def GME_SecondCalculateFreq(string):
		Result = {}
		Alt, Ref = [int(item) for item in string.split(',')]
		Result["AN"] = Alt + Ref
		Result["AC"] = Alt
		Result["AF"] = 0.0 if Result["AN"] == 0 else (Result["AC"] / Result["AN"])
		return Result

def GME_OMIM_Phenotype(string):
	shards = [item.strip() for item in string.split(';') if len(item) > 3]
	shards = [item[:-3].strip() for item in shards if item[-3:] in ["(1)", "(2)", "(3)", "(4)"]]
	shards = [{"MIM": None, "Desc": item.strip(' ,') } if not item[-6:].isdigit() else {"MIM": int(item[-6:]), "Desc": item[:-6].strip(' ,') } for item in shards]
	return shards

GME_Functions = {
	"str": lambda x: None if x in GME_Scheme["NoneValues"] else str(x),
	"int": lambda x: None if x in GME_Scheme["NoneValues"] else int(x),
	"int_comma": lambda x: None if x in GME_Scheme["NoneValues"] else [int(item) for item in x.split(',')],
	"str_comma": lambda x: None if x in GME_Scheme["NoneValues"] else [item for item in x.split(',')],
	"float": lambda x: None if x in GME_Scheme["NoneValues"] else float(x),
	"calculate_freq": lambda x: None if x in GME_Scheme["NoneValues"] else GME_CalculateFrequency(x),
	"omim_mim": lambda x: None if x in GME_Scheme["NoneValues"] else (int(x) if '[' not in x else json.loads(x.replace(' ', ', '))),
	"omim_disease": lambda x: None if x in GME_Scheme["NoneValues"] else GME_OMIM_Phenotype(x),
	"second_calculate_freq": lambda x: None if x in GME_Scheme["NoneValues"] else GME_SecondCalculateFreq(x)
}

def GME_ParseLine(Row):
	KeyList = list(GME_Scheme["InfoTags"].keys())
	NamedRow = { KeyList[index]: str(item) for index, item in enumerate(Row) }
	NewInfo = {}
	for Key, Value in NamedRow.items():
		ParsedValue = GME_Functions[GME_Scheme["InfoTags"][Key]["Func"]](Value)
		if ParsedValue is None: continue
		NewInfo[tuple(GME_Scheme["InfoTags"][Key]["Path"])] = ParsedValue
	Result = RecursiveDict(NewInfo)
	return Result

def GME_VariantMatch(Item, Row):
	return (Item[0] == Row[0]) and (int(Item[1]) == int(Row[1])) and (Item[2] == Row[2]) and (Item[3] == Row[3])

def GME_QueryFunction(Item, TabixObj):
	TabixOutput = TabixObj.query(Item[0], Item[1], Item[1])
	Result = []
	for Row in TabixOutput:
		if GME_VariantMatch(Item = Item, Row = Row):
			try:
				Result.append(GME_ParseLine(Row = Row))
			except Exception as Exc:
				raise ParseLineError(f"Row: {Row}; Exception Info: {Exc}")
	return None if not Result else Result
