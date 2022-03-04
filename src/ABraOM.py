from .SharedFunctions import * 
from .TabixFunctions import *

# ABRaOM 60+ SABE609 [exomes]

ABraOm_Scheme = LoadJSON(os.path.join(CurrentDir(), '..', 'schemes', 'ABRaOM_60+_SABE609_exomes_hg19_scheme.json'))

ABraOm_Functions = {
	"str": lambda x: None if x in ABraOm_Scheme["NoneValues"] else str(x),
	"str_comma": lambda x: None if x in ABraOm_Scheme["NoneValues"] else [item for item in x.split(',')],
	"int": int,
	"float": float
}

def ABraOm_ParseLine(Row):
	KeyList = list(ABraOm_Scheme["InfoTags"].keys())
	NamedRow = { KeyList[index]: str(item) for index, item in enumerate(Row) }
	NewInfo = {}
	for Key, Value in NamedRow.items():
		ParsedValue = ABraOm_Functions[ABraOm_Scheme["InfoTags"][Key]["Func"]](Value)
		if ParsedValue is None: continue
		NewInfo[tuple(ABraOm_Scheme["InfoTags"][Key]["Path"])] = ParsedValue
	Result = RecursiveDict(NewInfo)
	return Result

def ABraOm_VariantMatch(Item, Row):
	return (Item[0] == Row[0]) and (int(Item[1]) == int(Row[1])) and (Item[2] == Row[2]) and (Item[3] == Row[3])

def ABraOm_QueryFunction(Item, TabixObj):
	TabixOutput = TabixObj.query(Item[0], Item[1] - 1, Item[1])
	Result = []
	for Row in TabixOutput:
		if ABraOm_VariantMatch(Item = Item, Row = Row):
			try:
				Result.append(ABraOm_ParseLine(Row = Row))
			except Exception as Exc:
				raise ParseLineError(f"Row: {Row}; Exception Info: {Exc}")
	return None if not Result else Result
