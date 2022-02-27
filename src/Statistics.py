from src.SharedFunctions import *


## ------======| FASTQC |======------

def FastQC(InputFastQ, OutputHTML, Size = 0, Threads = multiprocessing.cpu_count()):
	logging.info(f'FastQ: "{InputFastQ}"; Report: "{OutputHTML}"')
	with tempfile.TemporaryDirectory() as TempDir:
		AnalyzeFilename = InputFastQ
		if Size != 0:
			logging.info(f'Subsample: {Size}')
			SampleFilename = os.path.join(TempDir, 'sample.fastq.gz')
			BashSubprocess(Name = f'FastQC.Sampling', Command = f'zcat -q "{InputFastQ}" | head -{Size * 4} | gzip -c > "{SampleFilename}"')
			AnalyzeFilename = SampleFilename
		BashSubprocess(Name = f'FastQC.Analysis', Command = f'fastqc -o "{TempDir}" -t {Threads} "{AnalyzeFilename}"')
		HTMLTemp = glob.glob(os.path.join(TempDir, '*.html'))
		if len(HTMLTemp) != 1:
			logging.error(f'Error processing file "{InputFastQ}"')
			raise RuntimeError
		BashSubprocess(Name = f'FastQC.Move', Command = f'cp "{HTMLTemp[0]}" "{OutputHTML}"') 


## ------======| FLAGSTAT |======------

def FlagStat(InputBAM, OutputJSON, Threads = multiprocessing.cpu_count()):
	BashSubprocess(Name = f'FlagStat.Statistics', Command = f'samtools flagstat -@ {Threads} -O json "{InputBAM}" > "{OutputJSON}"')

def LoadFlagStat(StatJSON): return LoadJSON(StatJSON)


## ------======| MARK DUPLICATES STATS |======------

def LoadMarkDuplicatesStat(StatTXT):
	Data = ''.join(open(StatTXT, 'rt').readlines())
	Stream = io.StringIO(Data.split('\n\n')[1])
	Table = pandas.read_csv(Stream, sep='\t', comment='#').set_index('LIBRARY')
	return Table.astype(object).transpose().to_dict()


## ------======| CUTADAPT STATS |======------

def CutadaptStat(StatTXT):
	Data = pandas.read_csv(StatTXT, sep='\t')
	return Data.transpose()[0].to_dict()


## ------======| COVERAGE & ENRICHMENT STATS |======------

def CoverageStats(Name, FinalBAM, StatsTXT, CaptureInfo, GenomeInfo):
	logging.info(f'BAM File: "{FinalBAM}"; Stats File: "{StatsTXT}"')
	with tempfile.TemporaryDirectory() as TempDir:
		CaptureTemp = os.path.join(TempDir, 'capture.csv')
		NotCaptureTemp = os.path.join(TempDir, 'not_capture.csv')
		BashSubprocess(Name = f'CoverageStats.Capture', Command = f'bedtools coverage -hist -sorted -g "{GenomeInfo["FAI"]}" -a "{CaptureInfo["CAP"]}" -b "{FinalBAM}" | grep -P "^all.*$" > "{CaptureTemp}"')
		BashSubprocess(Name = f'CoverageStats.NotCapture', Command = f'bedtools coverage -hist -sorted -g "{GenomeInfo["FAI"]}" -a "{CaptureInfo["CAP"]}" -b "{FinalBAM}" | grep -P "^all.*$" > "{NotCaptureTemp}"')
		CaptureData = pandas.read_csv(CaptureTemp, sep='\t', header=None, dtype={1: int, 4: float})[[1, 4]]
		NotCaptureData = pandas.read_csv(NotCaptureTemp, sep='\t', header=None, dtype={1: int, 4: float})[[1, 4]]
	Result = {
		"Name":                Name,
		"Capture DP>10 [%]":   CaptureData[CaptureData[1] >= 10][4].sum() * 100,
		"Capture Average":     CaptureData.apply(lambda x: x[1] * x[4], axis=1).sum(),
		"NotCapture Average":  NotCaptureData.apply(lambda x: x[1] * x[4], axis=1).sum(),
		"Enrichment Average":  None,
		"Capture DP0 [%]":     CaptureData.loc[0, 4] * 100,
		"NotCapture DP0 [%]":  NotCaptureData.loc[0, 4] * 100
	}
	try:
		Result["Enrichment Average"] = Result["Capture Average"] / Result["NotCapture Average"]
	except ZeroDivisionError:
		pass
	SaveJSON(Result, StatsTXT)

def LoadCoverageStats(StatJSON): return LoadJSON(StatJSON)
