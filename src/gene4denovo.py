from .SharedFunctions import * 
from .TabixFunctions import *

# gene4denovo

gene4denovo_Scheme = LoadJSON(os.path.join(CurrentDir(), '..', 'schemes', 'gene4denovo_v1.1_scheme.json'))

gene4denovo_Functions = {
	"str": lambda x: None if x in gene4denovo_Scheme["NullValues"] else str(x),
	"int": int
}

def gene4denovo_ParseLine(Row):
	KeyList = list(gene4denovo_Scheme["InfoTags"].keys())
	NamedRow = { KeyList[index]: str(item) for index, item in enumerate(Row) }
	NewInfo = {}
	for Key, Value in NamedRow.items():
		ParsedValue = gene4denovo_Functions[gene4denovo_Scheme["InfoTags"][Key]["Func"]](Value)
		if ParsedValue is None: continue
		NewInfo[tuple(gene4denovo_Scheme["InfoTags"][Key]["Path"])] = ParsedValue
	Result = RecursiveDict(NewInfo)
	return Result

def gene4denovo_VariantMatch(Item, Row):
	return (Item[0] == Row[0]) and (int(Item[1]) == int(Row[1])) and (Item[2] == Row[3]) and (Item[3] == Row[4])

def gene4denovo_QueryFunction(Item, TabixObj):
	TabixOutput = TabixObj.query(Item[0], Item[1], Item[1])
	Result = []
	for Row in TabixOutput:
		if gene4denovo_VariantMatch(Item = Item, Row = Row):
			try:
				Result.append(gene4denovo_ParseLine(Row = Row))
			except Exception as Exc:
				raise ParseLineError(f"Row: {Row}; Exception Info: {Exc}")
	return None if not Result else Result
