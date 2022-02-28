from .SharedFunctions import *


## ------======| TABIX |======------

def TabixThread(Chunk, FileName, PrepareFunction = None):
	ThreadNumber, Variants = Chunk
	StartTime = time.time()
	TabixObj = tabix.open(FileName)
	Result = []
	for index, item in Variants:
		try:
			Info = TabixObj.query(item[0], item[1], item[2])
		except tabix.TabixError:
			Info = None
		if PrepareFunction is not None: Info = PrepareFunction(Info)
		Result += [(index, Info)]
	logging.info(f'Name: TabixThread[{ThreadNumber}]; total time: {Timestamp(StartTime)}')
	return Result

def TabixSeek(List, FileName, PrepareFunction = None, Threads = multiprocessing.cpu_count()):
	ListArray = numpy.array_split(numpy.array(List, dtype=object), Threads)
	with Threading(f'Tabix.Seek', Threads = Threads):
		Result = sum(Pool.map(functools.partial(TabixThread, FileName = FileName, PrepareFunction = PrepareFunction), enumerate(ListArray)), [])
	return Result 
