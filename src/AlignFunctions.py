from .SharedFunctions import *


## ------======| MISC |======------

def MultipleTags(Tag, List, Quoted = True): return ' '.join([(f'{Tag} "{item}"' if Quoted else f'{Tag} {item}') for item in List])

def RGTag(Instrument, Lane, Platform, Barcode, Sample, Library):
	ID = f"D-{Instrument}.L-{Lane}"
	PL = Platform
	PU = f"D-{Instrument}.L-{Lane}.BC-{Barcode}"
	LB = f"LIB-{Sample}-{Library}"
	SM = Sample
	return f"@RG\\tID:{ID}\\tPL:{PL}\\tPU:{PU}\\tLB:{LB}\\tSM:{SM}"

## ------======| TRIM & CONVERT FASTQ |======------

def Solid2Illumina(InputFQ, OutputFQ, Threads = multiprocessing.cpu_count()):
	logging.info(f'Input FASTQ: {InputFQ}; Output FASTQ: {OutputFQ}')
	BashSubprocess(
		Name = f'Solid2Illumina.Convert',
		Command = f'cutadapt -j {Threads} -c --format=sra-fastq --bwa --action=none -o "{OutputFQ}" "{InputFQ}"'
	)
	
def Cutadapt(InputR1, InputR2, OutputR1, OutputR2, Adapter, ReportTXT, Threads = multiprocessing.cpu_count()):
	if InputR2 is None:
		logging.info(f'Mode: Single-end; Input FASTQ: "{InputR1}"; Output FASTQ: "{OutputR1}"; Adapter: {Adapter["Name"]}')
		Command = f'cutadapt --report minimal -j {Threads} -e 0.2 -m 8 -a {Adapter["R1"]} -o "{OutputR1}" "{InputR1}" > "{ReportTXT}"'
	else:
		logging.info(f'Mode: Single-end; Input FASTQ: "{InputR1}", "{InputR2}"; Output FASTQ: "{OutputR1}", "{OutputR2}"; Adapter: {Adapter["Name"]}')
		Command = f'cutadapt --report minimal -j {Threads} -e 0.2 -m 8 -a {Adapter["R1"]} -A {Adapter["R2"]} -o "{OutputR1}" -p "{OutputR2}" "{InputR1}" "{InputR2}" > "{ReportTXT}"'
	BashSubprocess(Name = f'Cutadapt.Trim', Command = Command)


## ------======| ALIGN & MERGE BAM |======------

def BWA(InputR1, InputR2, Reference, RGHeader, OutputBAM, Threads = multiprocessing.cpu_count()):
	if InputR2 is None:
		logging.info(f'Input FASTQ: "{InputR1}"; Output BAM: "{OutputBAM}"; Reference: "{Reference}"; RG Header: "{RGHeader}"')
		Command = ' '.join([
			f'set -o pipefail;',
			f'bwa mem -R "{RGHeader}" -t {Threads} -v 1 "{Reference}" "{InputR1}" |',
			f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["SortSam"]}" SortSam --VERBOSITY ERROR -SO queryname -I "/dev/stdin" -O "{OutputBAM}"'
		])
	else:
		logging.info(f'Input FASTQ: "{InputR1}", "{InputR2}"; Output BAM: "{OutputBAM}"; Reference: "{Reference}"; RG Header: "{RGHeader}"')
		Command = ' '.join([
			f'set -o pipefail;',
			f'bwa mem -R "{RGHeader}" -t {Threads} -v 1 "{Reference}" "{InputR1}" "{InputR2}" |',
			f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["SortSam"]}" SortSam --VERBOSITY ERROR -SO queryname -I "/dev/stdin" -O "{OutputBAM}"'
		])
	BashSubprocess(Name = f'BWA.Align&Sort', Command = Command)

def MergeBAMs(BAMs, OutputBAM, SortOrder):
	TaggedBAMs = MultipleTags(Tag = '-I', List = BAMs)
	BashSubprocess(
		Name = f'MergeBAMs.Merge',
		Command = f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["MergeSamFiles"]}" MergeSamFiles --VERBOSITY ERROR --USE_THREADING true -SO {SortOrder} {TaggedBAMs} -O "{OutputBAM}"'
	)


## ------======| MARK DUPLICATES |======------

def MarkDuplicates(InputBAM, OutputBAM, QuerySortedBAM, MetricsTXT, ActiveContigsTXT):
	logging.info(f'Input: "{InputBAM}"; Output: "{OutputBAM}"; Metrics: "{MetricsTXT}"')
	ActiveContigsCommand = ' '.join([
		f"samtools view -O SAM /dev/stdin |",
		f"awk -F '\\t' '{{ print $3 }}' - |",
		f"sort | uniq > \"{ActiveContigsTXT}\";"
	])
	BashSubprocess(
		Name = f'MarkDuplicates.Mark',
		Command = ' '.join([
			f'set -o pipefail;',
			f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["MarkDuplicates"]}" MarkDuplicates --REMOVE_DUPLICATES false --VERBOSITY ERROR --ADD_PG_TAG_TO_READS false --TAGGING_POLICY All --ASSUME_SORT_ORDER queryname -M "{MetricsTXT}" -I "{InputBAM}" -O "{QuerySortedBAM}"'
		])
	)
	BashSubprocess(
		Name = f'MarkDuplicates.RemoveAndSort',
		Command = ' '.join([
			f'samtools view -h -F 1024 "{QuerySortedBAM}" |',
			f'tee >( {ActiveContigsCommand} ) |',
			f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["SortSam"]}" SortSam --CREATE_INDEX true --VERBOSITY ERROR -SO coordinate -I "/dev/stdin" -O "{OutputBAM}"'
			])
		)

## ------======| BQSR |======------

def ContigBaseRecalibration(Contig, InputBAM, TempDir, dbSNP, Reference):
	FiltersComparison = [
		'MappedReadFilter',
		'MappingQualityAvailableReadFilter',
		'MappingQualityNotZeroReadFilter',
		'NotDuplicateReadFilter',
		'NotSecondaryAlignmentReadFilter',
		'PassesVendorQualityCheckReadFilter'
	]
	TaggedFilters = MultipleTags(Tag = '-RF', List = FiltersComparison, Quoted = False)
	# This is a bug above. BaseRecalibrator and ApplyBQSR filter reads differently, so ApplyBQSR f***s up every time it can't find RG.
	BQSRTable = os.path.join(TempDir, f'bqsr_table_{Contig}.tsv')
	OutputBAM = os.path.join(TempDir, f'output_{Contig}.bam')
	BashSubprocess(
		Name = f'ContigBQSR.MakeTable[{Contig}]',
		Command = f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["BaseRecalibrator"]}" BaseRecalibrator --verbosity ERROR -L {Contig} -I "{InputBAM}" --known-sites "{dbSNP}" -O "{BQSRTable}" -R "{Reference}"'
	)
	BashSubprocess(
		Name = f'ContigBQSR.Apply[{Contig}]',
		Command = f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["ApplyBQSR"]}" ApplyBQSR --verbosity ERROR {TaggedFilters} -OBI false -L "{Contig}" -bqsr "{BQSRTable}" -I "{InputBAM}" -O "{OutputBAM}"'
	)
	return OutputBAM

def BaseRecalibration(InputBAM, OutputBAM, dbSNP, Reference, ActiveContigs, Threads = multiprocessing.cpu_count()):
	logging.info(f'Input: "{InputBAM}"; Output: "{OutputBAM}"; Known sites: "{dbSNP}"; Reference: "{Reference}"')
	with tempfile.TemporaryDirectory() as TempDir:
		BashSubprocess(
			Name = f'BaseRecalibration.PreIndex',
			Command = f'{GATK_PATH} BuildBamIndex -I "{InputBAM}"'
		)
		with Threading('ContigBQSR', Threads) as pool:
			Shards = pool.map(functools.partial(ContigBaseRecalibration, InputBAM = InputBAM, TempDir = TempDir, dbSNP = dbSNP, Reference = Reference), ActiveContigs)
			TaggedShards = MultipleTags(Tag = '-I', List = Shards)
		BashSubprocess(
			Name = f'BaseRecalibration.Merge',
			Command = f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["MergeSamFiles"]}" MergeSamFiles --VERBOSITY ERROR --CREATE_MD5_FILE true --CREATE_INDEX true --USE_THREADING true -SO coordinate {TaggedShards} -O "{OutputBAM}"'
		)


## ------======| VARIANT CALLING |======------

def ContigHaplotypeCalling(Contig, InputBAM, TempDir, Reference):
	OutputGVCF = os.path.join(TempDir, f'output_{Contig}.g.vcf.gz')
	OutputVCF = os.path.join(TempDir, f'output_{Contig}.vcf.gz')
	BashSubprocess(
		Name = f'ContigCalling.Calling[{Contig}]',
		Command = f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["HaplotypeCaller"]}" HaplotypeCaller -OVI true --verbosity ERROR --native-pair-hmm-threads 2 --dont-use-soft-clipped-bases true -ERC GVCF -L "{Contig}" -I "{InputBAM}" -O "{OutputGVCF}" -R "{Reference}"'
	)
	BashSubprocess(
		Name = f'ContigCalling.Genotype[{Contig}]',
		Command = f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["GenotypeGVCFs"]}" GenotypeGVCFs --verbosity ERROR -L "{Contig}" -V "{OutputGVCF}" -O "{OutputVCF}" -R "{Reference}"'
	)
	return { 'VCF': OutputVCF, 'gVCF': OutputGVCF }

def HaplotypeCalling(InputBAM, OutputVCF, OutputGVCF, Reference, ActiveContigs, Threads = multiprocessing.cpu_count()):
	logging.info(f'Input BAM: "{InputBAM}"; Output VCF: "{OutputVCF}"; Output gVCF: "{OutputGVCF}"; Reference: "{Reference}"')
	with tempfile.TemporaryDirectory() as TempDir:
		with Threading('ContigCalling', Threads) as pool:
			Shards = pool.map(functools.partial(ContigHaplotypeCalling, InputBAM = InputBAM, TempDir = TempDir, Reference = Reference), ActiveContigs)
			TaggedShardsVCF = MultipleTags(Tag = '-I', List = [item['VCF'] for item in Shards])
			TaggedShardsGVCF = MultipleTags(Tag = '-I', List = [item['gVCF'] for item in Shards])
		BashSubprocess(
			Name = f'HaplotypeCalling.MergeVCF',
			Command = f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["MergeVcfs"]}" MergeVcfs --CREATE_INDEX true --VERBOSITY ERROR --CREATE_MD5_FILE true {TaggedShardsVCF} -O "{OutputVCF}"'
		)
		BashSubprocess(
			Name = f'HaplotypeCalling.MergeGVCF',
			Command = f'{GATK_PATH} --java-options "{CONFIG_JAVA_OPTIONS["MergeVcfs"]}" MergeVcfs --CREATE_INDEX true --VERBOSITY ERROR --CREATE_MD5_FILE true {TaggedShardsGVCF} -O "{OutputGVCF}"'
		)
