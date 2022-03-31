from .SharedFunctions import *


## ------======| FASTQC |======------

def FastQC(**kwargs):
	logging.info('*')
	for Key, Value in kwargs.items(): logging.info(RenderParameters(Key, Value))
	N = types.SimpleNamespace(**kwargs)
	with tempfile.TemporaryDirectory() as TempDir:
		AnalyzeFilename = N.Input_FastQ
		if 'Subsample_Size' in kwargs:
			SampleFilename = os.path.join(TempDir, 'sample.fastq.gz')
			CommandSampling = f'zcat -q "{N.Input_FastQ}" | head -{N.Subsample_Size * 4} | gzip -c > "{SampleFilename}"'
			BashSubprocess(Name = f'FastQC.Sampling', Command = CommandSampling)
			AnalyzeFilename = SampleFilename
		CommandAnalysis = f'fastqc -o "{TempDir}" -t {N.Threads} "{AnalyzeFilename}"'
		BashSubprocess(Name = f'FastQC.Analysis', Command = CommandAnalysis)
		HTMLTemp = glob.glob(os.path.join(TempDir, '*.html'))
		if len(HTMLTemp) != 1:
			logging.error(f'Error processing file "{N.Input_FastQ}"')
			raise RuntimeError
		CommandMove = f'mv "{HTMLTemp[0]}" "{N.Output_HTML}"'
		BashSubprocess(Name = f'FastQC.Move', Command = CommandMove) 


## ------======| FLAGSTAT |======------

def FlagStat(**kwargs):
	logging.info('*')
	for Key, Value in kwargs.items(): logging.info(RenderParameters(Key, Value))
	N = types.SimpleNamespace(**kwargs)
	Command = f'samtools flagstat -@ {N.Threads} -O json "{N.Input_BAM}" > "{N.Samtools_Flagstats}"'
	BashSubprocess(Name = f'FlagStat', Command = Command)


## ------======| MARK DUPLICATES STATS |======------

def LoadMarkDuplicatesStat(StatTXT):
	logging.info('*')
	Data = ''.join(open(StatTXT, 'rt').readlines())
	Stream = io.StringIO(Data.split('\n\n')[1])
	Table = pandas.read_csv(Stream, sep='\t', comment='#').set_index('LIBRARY')
	return Table.astype(object).transpose().to_dict()


## ------======| CUTADAPT STATS |======------

def CutadaptStat(StatTXT):
	logging.info('*')
	Data = pandas.read_csv(StatTXT, sep='\t')
	return Data.transpose()[0].to_dict()


## ------======| COVERAGE & ENRICHMENT STATS |======------

def CoverageStats(**kwargs):
	logging.info('*')
	for Key, Value in kwargs.items(): logging.info(RenderParameters(Key, Value))
	N = types.SimpleNamespace(**kwargs)
	LoadBedtoolsOutput = lambda FN: pandas.read_csv(FN, sep='\t', header=None, dtype={1: int, 4: float})[[1, 4]]
	with tempfile.TemporaryDirectory() as TempDir:
		CaptureTemp = os.path.join(TempDir, 'capture.csv')
		NotCaptureTemp = os.path.join(TempDir, 'not_capture.csv')
		CommandCaptureCoverage = f'bedtools coverage -hist -sorted -g "{N.Reference_FAI}" -a "{N.Capture_BED}" -b "{N.Input_BAM}" | grep -P "^all.*$" > "{CaptureTemp}"'
		CommandNotCaptureCoverage = f'bedtools coverage -hist -sorted -g "{N.Reference_FAI}" -a "{N.NotCapture_BED}" -b "{N.Input_BAM}" | grep -P "^all.*$" > "{NotCaptureTemp}"'
		BashSubprocess(Name = f'CoverageStats.Capture', Command = CommandCaptureCoverage)
		BashSubprocess(Name = f'CoverageStats.NotCapture', Command = CommandNotCaptureCoverage)
		CaptureData = LoadBedtoolsOutput(CaptureTemp)
		NotCaptureData = LoadBedtoolsOutput(NotCaptureTemp)
	Result = {
		'Capture DP>10 [%]':   CaptureData[CaptureData[1] >= 10][4].sum() * 100,
		'Capture Average':     CaptureData.apply(lambda x: x[1] * x[4], axis=1).sum(),
		'NotCapture Average':  NotCaptureData.apply(lambda x: x[1] * x[4], axis=1).sum(),
		'Enrichment Average':  None,
		'Capture DP0 [%]':     CaptureData.loc[0, 4] * 100,
		'NotCapture DP0 [%]':  NotCaptureData.loc[0, 4] * 100
	}
	try:
		Result['Enrichment Average'] = Result['Capture Average'] / Result['NotCapture Average']
	except ZeroDivisionError:
		pass
	SaveJSON(Result, N.Coverage_Stats)
