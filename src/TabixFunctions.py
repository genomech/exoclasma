from .SharedFunctions import *


## ------======| TABIX |======------

def TabixThread(Chunk, FileName, QueryFunction):
	ThreadNumber, Variants = Chunk
	StartTime = time.time()
	TabixObj = tabix.open(FileName)
	Result = []
	for index, item in Variants:
		Info = QueryFunction(item, TabixObj)
		Result += [(index, Info)]
	logging.info(f'Name: TabixThread[{ThreadNumber}]; total time: {Timestamp(StartTime)}')
	return Result

def TabixSeek(List, FileName, QueryFunction, Threads = multiprocessing.cpu_count()):
	ListArray = numpy.array_split(numpy.array(List, dtype=object), Threads)
	with Threading(f'Tabix.Seek', Threads = Threads) as Pool:
		Result = sum(Pool.map(functools.partial(TabixThread, FileName = FileName, QueryFunction = QueryFunction), enumerate(ListArray)), [])
	return Result 

## ------======| VCF ROW |======------

def ParseVcfRow(Row):
	try:
		Qual = int(Row[5])
	except ValueError:
		try:
			Qual = float(Row[5])
		except ValueError:
			Qual = None
	try:
		Format = Row[8]
	except IndexError:
		Format = None
	else:
		if Format == '.': Format = None
	Result = {
		"CHROM": Row[0],
		"POS": int(Row[1]),
		"ID": None if (Row[2] == '.') else Row[2],
		"REF": Row[3],
		"ALT": Row[4].split(','),
		"QUAL": Qual,
		"FILTER": None if (Row[6] == '.') else ([] if (Row[6] == 'PASS') else Row[6].split(';')),
		"INFO": Row[7]
	}
	Result["INFO"] = [i.split('=') for i in Result["INFO"].split(';')]
	Result["INFO"] = dict([(i[0], True) if (len(i) == 1) else ((i[0], i[1]) if (len(i) == 2) else None) for i in Result["INFO"]])
	return Result
	
