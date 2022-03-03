from .SharedFunctions import * 
from .TabixFunctions import *

# dbNSFP v4.1a

dbNSFP_Scheme = LoadJSON(os.path.join(CurrentDir(), '..', 'schemes', 'dbNSFP4.1a_scheme.json'))

FloatNone = lambda x: None if x in dbNSFP_Scheme['NoneValues'] else float(x)
IntNone = lambda x: None if x in dbNSFP_Scheme['NoneValues'] else int(x)
StrNone = lambda x: None if x in dbNSFP_Scheme['NoneValues'] else str(x)

dbNSFP_Functions = {
	'clinvar_source': lambda x: None if x in dbNSFP_Scheme['NoneValues'] else [{'DB': item.split(':', maxsplit=1)[0].replace('_', ' '), 'ID': item.split(':', maxsplit=1)[1]} for item in x.split('|')],
	'float': FloatNone,
	'float_semicolon': lambda x: None if x in dbNSFP_Scheme['NoneValues'] else [FloatNone(item) for item in x.split(';')],
	'float_semicolon_withoutlast': lambda x: None if x in dbNSFP_Scheme['NoneValues'] else [FloatNone(item) for item in x.split(';')[:-1]],
	'int': IntNone,
	'int_semicolon': lambda x: None if x in dbNSFP_Scheme['NoneValues'] else [IntNone(item) for item in x.split(';')], 
	'int_vline': lambda x: None if x in dbNSFP_Scheme['NoneValues'] else [IntNone(item) for item in x.split('|')],
	'siphy_pi': lambda x: None if x in dbNSFP_Scheme['NoneValues'] else { "ACGT"[index]: float(item) for index, item in enumerate(x.split(':'))},
	'str': StrNone,
	'str_semicolon': lambda x: None if x in dbNSFP_Scheme['NoneValues'] else [StrNone(item) for item in x.split(';')],
	'str_semicolon_vline': lambda x: None if x in dbNSFP_Scheme['NoneValues'] else [(None if item in dbNSFP_Scheme['NoneValues'] else [StrNone(item2) for item2 in item.split('|')]) for item in x.split(';')],
	'str_semicolon_withoutlast': lambda x: None if x in dbNSFP_Scheme['NoneValues'] else [StrNone(item) for item in x.split(';')[:-1]],
	'str_underline_purge': lambda x: StrNone(x.replace('_', ' ')),
	'str_vline': lambda x: None if x in dbNSFP_Scheme['NoneValues'] else [StrNone(item) for item in x.split('|')],
	'str_vline_underline_purge': lambda x: None if x in dbNSFP_Scheme['NoneValues'] else [StrNone(item.replace('_', ' ')) for item in x.split('|')],
	'top5features': lambda x: None if x in dbNSFP_Scheme['NoneValues'] else [{ ("Desc", "P")[index2]: item2 for index2, item2 in enumerate(item[:-1].split(' (P = ')) } for item in x.split('; ')]
}

def dbNSFP_ParseLine(Row):
	KeyList = list(dbNSFP_Scheme["InfoTags"].keys())
	NamedRow = { KeyList[index]: item for index, item in enumerate(Row) }
	NewInfo = {}
	for Key, Value in NamedRow.items():
		ParsedValue = dbNSFP_Functions[dbNSFP_Scheme["InfoTags"][Key]["Func"]](Value)
		if ParsedValue is None: continue
		if type(ParsedValue) == list:
			if all([item is None for item in ParsedValue]): continue
		NewInfo[tuple(dbNSFP_Scheme["InfoTags"][Key]["Path"])] = ParsedValue
	Result = RecursiveDict(NewInfo)
	return Result

def dbNSFP_VariantMatch(Item, Row):
	return (Item[0] == Row[7]) and (int(Item[1]) == int(Row[8])) and (Item[2] == Row[2]) and (Item[3] == Row[3]) # FIXME hg19 only

def dbNSFP_QueryFunction(Item, TabixObj):
	TabixOutput = TabixObj.query(Item[0], Item[1] - 1, Item[1])
	Result = [ dbNSFP_ParseLine(Row = Row) for Row in TabixOutput if dbNSFP_VariantMatch(Item = Item, Row = Row) ]
	return None if not Result else Result
