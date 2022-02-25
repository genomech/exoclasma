from src.SharedFunctions import *


## ------======| MISC |======------

def MultipleTags(Tag, List, Quoted = True): return ' '.join([(f'{Tag} "{item}"' if Quoted else f'{Tag} {item}') for item in List])

def GetActiveContigs(InputBAM, Threads = multiprocessing.cpu_count()): 
	BashSubprocess(Name = f'GetActiveContigs', Command = f'samtools view -@ {Threads} "{InputBAM}" | awk -F \'\\t\' \'{{ print $3 }}\' - | uniq').decode('utf-8')[:-1].split('\n')


## ------======| TRIM & CONVERT FASTQ |======------

def Solid2Illumina(InputFQ, OutputFQ, Threads = multiprocessing.cpu_count()):
	logging.info(f'Input FASTQ: {InputFQ}; Output FASTQ: {OutputFQ}')
	BashSubprocess(Name = f'Solid2Illumina.Convert', Command = f'cutadapt -j {Threads} -c --format=sra-fastq --bwa --action=none -o "{OutputFQ}" "{InputFQ}"')
	
def Cutadapt(InputR1, InputR2, OutputR1, OutputR2, Adapter, ReportTXT, Threads = multiprocessing.cpu_count()):
	if InputR2 is None:
		logging.info(f'Mode: Single-end; Input FASTQ: "{InputR1}"; Output FASTQ: "{OutputR1}"; Adapter: {Adapter["Name"]}')
		Command = f'cutadapt -j {Threads} -e 0.2 -m 8 -a {Adapter["R1"]} -o "{OutputR1}" "{InputR1}" > "{ReportTXT}"'
	else:
		logging.info(f'Mode: Single-end; Input FASTQ: "{InputR1}", "{InputR2}"; Output FASTQ: "{OutputR1}", "{OutputR2}"; Adapter: {Adapter["Name"]}')
		Command = f'cutadapt -j {Threads} -e 0.2 -m 8 -a {Adapter["R1"]} -A {Adapter["R2"]} -o "{OutputR1}" -p "{OutputR2}" "{InputR1}" "{InputR2}" > "{ReportTXT}"'
	BashSubprocess(Name = f'Cutadapt.Trim', Command = Command)


## ------======| ALIGN & MERGE BAM |======------

def BWA(InputR1, InputR2, Reference, RGHeader, OutputBAM, Threads = multiprocessing.cpu_count()):
	if InputR2 is None:
		logging.info(f'Input FASTQ: "{InputR1}"; Output BAM: "{OutputBAM}"; Reference: "{Reference}"; RG Header: "{RGHeader}"')
		Command = f'set -o pipefail; bwa mem -R "{RGHeader}" -t {Threads} -v 1 "{Reference}" "{InputR1}" | {GATK_PATH} SortSam --VERBOSITY ERROR -SO queryname -I "/dev/stdin" -O "{OutputBAM}"'
	else:
		logging.info(f'Input FASTQ: "{InputR1}", "{InputR2}"; Output BAM: "{OutputBAM}"; Reference: "{Reference}"; RG Header: "{RGHeader}"')
		Command = f'bwa mem -R "{RGHeader}" -t {Threads} -v 1 "{Reference}" "{InputR1}" "{InputR2}" | {GATK_PATH} SortSam --VERBOSITY ERROR -SO queryname -I "/dev/stdin" -O "{OutputBAM}"'
	BashSubprocess(Name = f'BWA.Align&Sort', Command = Command)

def MergeBAMs(BAMs, OutputBAM, SortOrder):
	TaggedBAMs = MultipleTags(Tag = '-I', List = BAMs)
	BashSubprocess(Name = f'MergeBAMs.Merge', Command = f'{GATK_PATH} MergeSamFiles --USE_THREADING true -SO {SortOrder} {TaggedBAMs} -O "{OutputBAM}"')


## ------======| MARK DUPLICATES |======------

def MarkDuplicates(InputBAM, OutputBAM, MetricsTXT):
	logging.info(f'Input: "{InputBAM}"; Output: "{OutputBAM}"; Metrics: "{MetricsTXT}"')
	BashSubprocess(Name = f'MarkDuplicates.RemoveAndSort', Command = f'set -o pipefail; {GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["MarkDuplicates"]}" MarkDuplicates --REMOVE_DUPLICATES true --VERBOSITY ERROR --ASSUME_SORT_ORDER queryname -M "{MetricsTXT}" -I "{InputBAM}" -O "/dev/stdout" | {GATK_PATH} SortSam --VERBOSITY ERROR -SO coordinate -I "/dev/stdin" -O "{OutputBAM}"')
	BashSubprocess(Name = f'MarkDuplicates.Index', Command = f'{GATK_PATH} BuildBamIndex -I "{OutputBAM}"')


## ------======| BQSR |======------

def ContigBaseRecalibration(Contig, InputBAM, TempDir, dbSNP, Reference):
	FiltersComparison = ['MappedReadFilter', 'MappingQualityAvailableReadFilter', 'MappingQualityNotZeroReadFilter', 'NotDuplicateReadFilter', 'NotSecondaryAlignmentReadFilter', 'PassesVendorQualityCheckReadFilter']
	TaggedFilters = MultipleTags(Tag = '-RF', List = FiltersComparison, Quoted = False)
	# This is a bug above. BaseRecalibrator and ApplyBQSR filter reads differently, so ApplyBQSR f***s up every time it can't find RG.
	BQSRTable = os.path.join(TempDir, f'bqsr_table_{Contig}.tsv')
	OutputBAM = os.path.join(TempDir, f'output_{Contig}.bam')
	BashSubprocess(Name = f'ContigBQSR.MakeTable[{Contig}]', Command = f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["BaseRecalibrator"]}" BaseRecalibrator --tmp-dir "{TempDir}" -L {Contig} -I "{InputBAM}" --known-sites "{dbSNP}" -O "{BQSRTable}" -R "{Reference}"')
	BashSubprocess(Name = f'ContigBQSR.Apply[{Contig}]', Command = f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["ApplyBQSR"]}" ApplyBQSR {TaggedFilters} --tmp-dir "{TempDir}" -OBI false -L "{Contig}" -bqsr "{BQSRTable}" -I "{InputBAM}" -O "{OutputBAM}"')
	return OutputBAM

def BaseRecalibration(InputBAM, OutputBAM, dbSNP, Reference, ActiveContigs, Threads = multiprocessing.cpu_count()):
	logging.info(f'Input: "{InputBAM}"; Output: "{OutputBAM}"; Known sites: "{dbSNP}"; Reference: "{Reference}"')
	with tempfile.TemporaryDirectory() as TempDir:
		BashSubprocess(Name = f'BaseRecalibration.PreIndex', Command = f'{GATK_PATH} BuildBamIndex -I "{InputBAM}"')
		with Threading('ContigBQSR', Threads) as pool:
			Shards = pool.map(functools.partial(ContigBaseRecalibration, InputBAM = InputBAM, TempDir = TempDir, dbSNP = dbSNP, Reference = Reference), ActiveContigs)
			TaggedShards = MultipleTags(Tag = '-I', List = Shards)
		BashSubprocess(Name = f'BaseRecalibration.Merge', Command = f'{GATK_PATH} MergeSamFiles --USE_THREADING true -SO coordinate {TaggedShards} -O "{OutputBAM}"')
		BashSubprocess(Name = f'BaseRecalibration.PostIndex', Command = f'{GATK_PATH} BuildBamIndex -I "{OutputBAM}"')


## ------======| VARIANT CALLING |======------

def ContigHaplotypeCalling(Contig, InputBAM, TempDir, Reference):
	OutputVCF = os.path.join(TempDir, f'output_{Contig}.vcf')
	BashSubprocess(Name = f'ContigCalling.Calling[{Contig}]', Command = f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["HaplotypeCalling"]}" HaplotypeCaller --native-pair-hmm-threads 2 -OVI false --dont-use-soft-clipped-bases true -L "{Contig}" -I "{InputBAM}" -O "{OutputVCF}" -R "{Reference}"')
	return OutputVCF

def HaplotypeCalling(InputBAM, OutputVCF, Reference, ActiveContigs, Threads = multiprocessing.cpu_count()):
	logging.info(f'Input BAM: "{InputBAM}"; Output VCF: "{OutputVCF}"; Reference: "{Reference}"')
	with tempfile.TemporaryDirectory() as TempDir:
		with Threading('ContigCalling', Threads) as pool:
			Shards = pool.map(functools.partial(ContigHaplotypeCalling, InputBAM = InputBAM, TempDir = TempDir, Reference = Reference), ActiveContigs)
			TaggedShards = MultipleTags(Tag = '-I', List = Shards)
		BashSubprocess(Name = f'HaplotypeCalling.Merge', Command = f'{GATK_PATH} MergeVcfs {TaggedShards} -O "{OutputVCF}"')
