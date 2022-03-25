from .SharedFunctions import *


## ------======| TABIX |======------

class ParseLineError(Exception): pass

def TabixThread(Chunk, FileName, QueryFunction):
	ThreadNumber, Variants = Chunk
	StartTime = time.time()
	TabixObj = tabix.open(FileName)
	Result = []
	for index, item in Variants:
		Info = QueryFunction(item, TabixObj)
		Result.append((index, Info))
	logging.info(f'Name: TabixThread[{ThreadNumber}]; total time: {Timestamp(StartTime)}')
	return Result

def TabixSeek(List, FileName, QueryFunction, Threads = multiprocessing.cpu_count()):
	ListArray = numpy.array_split(numpy.array(List, dtype=object), Threads)
	with Threading(f'Tabix.Seek', Threads = Threads) as Pool:
		Result = sum(Pool.map(functools.partial(TabixThread, FileName = FileName, QueryFunction = QueryFunction), enumerate(ListArray)), [])
	return Result 

## ------======| VCF ROW |======------

def ExtractSamples(Vcf):
	Stream = OpenAnyway(Vcf, 'rt')
	while 1:
		Row = next(Stream)
		if Row[:6] == "#CHROM": break
	return {index: item for index, item in enumerate(Row[:-1].split('\t')[9:])}

def ParseVcfRow(Row, Samples = None):
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
	INFO = [i.split('=') for i in Result["INFO"].split(';')]
	Result["INFO"] = {}
	Tags = [item[0] for item in INFO]
	TagsCount = {item[0]: 0 for item in INFO}
	for item in INFO:
		if Tags.count(item[0]) > 1:
			tag = tuple([item[0], TagsCount[item[0]]])
			TagsCount[item[0]] += 1
		else: tag = item[0]
		Result["INFO"][tag] = True if len(item) == 1 else item[1]
	if len(Row) > 8:
		assert Samples is not None, f"No samples!"
		Format = Row[8].split(':')
		Result["SAMPLES"] = { Samples[index]: { Format[i]: value for i, value in enumerate(item.split(':'))} for index, item in enumerate(Row[9:]) }
	return Result

def VcfVariantMatch(Item, Row):
	return (Item[0] == Row[0]) and (int(Item[1]) == int(Row[1])) and (Item[2] == Row[3]) and (Item[3] in Row[4].split(','))
