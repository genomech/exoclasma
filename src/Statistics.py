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

def FlagStat(InputBAM, OutputJSON, Threads = multiprocessing.cpu_count()):
	BashSubprocess(Name = f'FlagStat.Statistics', Command = f'samtools flagstat -@ {Threads} -O json "{InputBAM}" > "{OutputJSON}"')
