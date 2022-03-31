from .SharedFunctions import *


## ------======| MISC |======------

def MultipleTags(Tag, List, Quoted = True): return ' '.join([(f'{Tag} "{item}"' if Quoted else f'{Tag} {item}') for item in List])

def RGTag(**kwargs):
	logging.info('* entry point *')
	N = types.SimpleNamespace(**kwargs)
	ID = f'D-{N.Instrument}.L-{N.Lane}'
	PL = N.Platform
	PU = f'D-{N.Instrument}.L-{N.Lane}.BC-{N.Barcode}'
	LB = f'LIB-{N.Sample}-{N.Library}'
	SM = N.Sample
	return f'@RG\\tID:{ID}\\tPL:{PL}\\tPU:{PU}\\tLB:{LB}\\tSM:{SM}'

## ------======| TRIM & CONVERT FASTQ |======------

def SolidToIllumina(**kwargs):
	logging.info('* entry point *')
	N = RenderParameters(kwargs)
	Command = f'cutadapt -j {N.Threads} -c --format=sra-fastq --bwa --action=none -o "{N.Output_FastQ}" "{N.Input_FastQ}"'
	BashSubprocess(Name = 'SolidToIllumina', Command = Command)
	
def Cutadapt(**kwargs):
	logging.info('* entry point *')
	N = RenderParameters(kwargs)
	if N.Mode == 'Single-end':
		Command = f'cutadapt --report minimal -j {N.Threads} -e 0.2 -m 8 -a {CONFIG_ADAPTERS[N.Adapter]["R1"]} -o "{N.Output_FastQ}" "{N.Input_FastQ}" > "{N.Cutadapt_Report}"'
	if N.Mode == 'Paired-end':
		Command = f'cutadapt --report minimal -j {N.Threads} -e 0.2 -m 8 -a {CONFIG_ADAPTERS[N.Adapter]["R1"]} -A {CONFIG_ADAPTERS[N.Adapter]["R2"]} -o "{N.Output_FastQ_R1}" -p "{N.Output_FastQ_R2}" "{N.Input_FastQ_R1}" "{N.Input_FastQ_R2}" > "{N.Cutadapt_Report}"'
	BashSubprocess(Name = 'Cutadapt', Command = Command)


## ------======| ALIGN & MERGE BAM |======------

def BWA(**kwargs):
	logging.info('* entry point *')
	N = RenderParameters(kwargs)
	if N.Mode == 'Single-end':
		Command = f'set -o pipefail; bwa mem -R "{N.RG_Header}" -t {N.Threads} -v 1 "{N.Reference_Fasta}" "{N.Input_FastQ}" | "{GATK_PATH}" --java-options "{CONFIG_JAVA_OPTIONS["SortSam"]}" SortSam --VERBOSITY ERROR -SO queryname -I "/dev/stdin" -O "{N.Output_BAM}"'
	if N.Mode == 'Paired-end':
		Command = f'set -o pipefail; bwa mem -R "{N.RG_Header}" -t {N.Threads} -v 1 "{N.Reference_Fasta}" "{N.Input_FastQ_R1}" "{N.Input_FastQ_R2}" | "{GATK_PATH}" --java-options "{CONFIG_JAVA_OPTIONS["SortSam"]}" SortSam --VERBOSITY ERROR -SO queryname -I "/dev/stdin" -O "{N.Output_BAM}"'
	BashSubprocess(Name = 'BWA', Command = Command)

def MergeSamFiles(**kwargs):
	logging.info('* entry point *')
	N = RenderParameters(kwargs)
	TaggedBAMs = MultipleTags(Tag = '-I', List = N.BAM_List)
	Command = f'"{GATK_PATH}" --java-options "{CONFIG_JAVA_OPTIONS["MergeSamFiles"]}" MergeSamFiles --VERBOSITY ERROR --USE_THREADING true -SO {N.Sort_Order} {TaggedBAMs} -O "{N.Output_BAM}"'
	BashSubprocess(Name = 'MergeSamFiles', Command = Command)


## ------======| MARK DUPLICATES |======------

def MarkDuplicates(**kwargs):
	logging.info('* entry point *')
	N = RenderParameters(kwargs)
	ActiveContigsCommand = f'samtools view -O SAM "/dev/stdin" | awk -F \'\\t\' \'{{ print $3 }}\' - | sort | uniq > "{N.Active_Contigs_File}";'
	CommandMarkDuplicates = f'set -o pipefail; "{GATK_PATH}" --java-options "{CONFIG_JAVA_OPTIONS["MarkDuplicates"]}" MarkDuplicates --REMOVE_DUPLICATES false --VERBOSITY ERROR --ADD_PG_TAG_TO_READS false --TAGGING_POLICY All --ASSUME_SORT_ORDER queryname -M "{N.MarkDuplicates_Metrics}" -I "{N.Input_BAM}" -O "{N.Query_Sorted_BAM}"'
	CommandRemoveAndSort = f'samtools view -h -F 1024 "{N.Query_Sorted_BAM}" | tee >( {ActiveContigsCommand} ) | "{GATK_PATH}" --java-options "{CONFIG_JAVA_OPTIONS["SortSam"]}" SortSam --CREATE_INDEX true --VERBOSITY ERROR -SO coordinate -I "/dev/stdin" -O "{N.Output_BAM}"'
	BashSubprocess(Name = 'MarkDuplicates.MarkDups', Command = CommandMarkDuplicates)
	BashSubprocess(Name = 'MarkDuplicates.RemoveDupsAndSortSam', Command = CommandRemoveAndSort)

## ------======| BQSR |======------

def ContigBaseRecalibration(Contig, **kwargs):
	N = types.SimpleNamespace(**kwargs)
	FiltersComparison = [
		'MappedReadFilter',
		'MappingQualityAvailableReadFilter',
		'MappingQualityNotZeroReadFilter',
		'NotDuplicateReadFilter',
		'NotSecondaryAlignmentReadFilter',
		'PassesVendorQualityCheckReadFilter'
	]
	TaggedFilters = MultipleTags(Tag = '-RF', List = FiltersComparison, Quoted = False)
	BQSRTable = os.path.join(N.Temp_Directory, f'bqsr_table_{Contig}.tsv')
	OutputBAM = os.path.join(N.Temp_Directory, f'output_{Contig}.bam')
	CommandBaseRecalibrator = f'"{GATK_PATH}" --java-options "{CONFIG_JAVA_OPTIONS["BaseRecalibrator"]}" BaseRecalibrator --verbosity ERROR -L {Contig} -I "{N.Input_BAM}" --known-sites "{N.dbSNP_Known_Sites}" -O "{BQSRTable}" -R "{N.Reference_Fasta}"'
	CommandApplyBQSR = f'"{GATK_PATH}" --java-options "{CONFIG_JAVA_OPTIONS["ApplyBQSR"]}" ApplyBQSR --verbosity ERROR {TaggedFilters} -OBI false -L "{Contig}" -bqsr "{BQSRTable}" -I "{N.Input_BAM}" -O "{OutputBAM}"'
	BashSubprocess(Name = f'BaseRecalibration.Contig.{Contig}.BaseRecalibrator', Command = CommandBaseRecalibrator)
	BashSubprocess(Name = f'BaseRecalibration.Contig.{Contig}.ApplyBQSR', Command = CommandApplyBQSR)
	return OutputBAM

def BaseRecalibration(**kwargs):
	logging.info('* entry point *')
	N = RenderParameters(kwargs)
	with tempfile.TemporaryDirectory() as TempDir:
		CommandBuildBamIndex = f'"{GATK_PATH}" BuildBamIndex -I "{N.Input_BAM}"'
		BashSubprocess(Name = 'BaseRecalibration.BuildBamIndex', Command = CommandBuildBamIndex)
		with Threading('BaseRecalibration.Contig', N.Threads) as pool:
			Shards = pool.map(functools.partial(ContigBaseRecalibration, Input_BAM = N.Input_BAM, Temp_Directory = TempDir, dbSNP_Known_Sites = N.dbSNP_Known_Sites, Reference_Fasta = N.Reference_Fasta), N.Active_Contigs)
			TaggedShards = MultipleTags(Tag = '-I', List = Shards)
		CommandMergeSamFiles = f'"{GATK_PATH}" --java-options "{CONFIG_JAVA_OPTIONS["MergeSamFiles"]}" MergeSamFiles --VERBOSITY ERROR --CREATE_MD5_FILE true --CREATE_INDEX true --USE_THREADING true -SO coordinate {TaggedShards} -O "{N.Output_BAM}"'
		BashSubprocess(Name = 'BaseRecalibration.MergeSamFiles', Command = CommandMergeSamFiles)


## ------======| VARIANT CALLING |======------

def ContigHaplotypeCalling(Contig, **kwargs):
	N = types.SimpleNamespace(**kwargs)
	OutputGVCF = os.path.join(N.Temp_Directory, f'output_{Contig}.g.vcf.gz')
	OutputVCF = os.path.join(N.Temp_Directory, f'output_{Contig}.vcf.gz')
	CommandHaplotypeCaller = f'"{GATK_PATH}" --java-options "{CONFIG_JAVA_OPTIONS["HaplotypeCaller"]}" HaplotypeCaller -OVI true --verbosity ERROR --native-pair-hmm-threads 2 --dont-use-soft-clipped-bases true -ERC GVCF -L "{Contig}" -I "{N.Input_BAM}" -O "{OutputGVCF}" -R "{N.Reference_Fasta}"'
	CommandGenotypeGVCFs = f'"{GATK_PATH}" --java-options "{CONFIG_JAVA_OPTIONS["GenotypeGVCFs"]}" GenotypeGVCFs --verbosity ERROR -L "{Contig}" -V "{OutputGVCF}" -O "{OutputVCF}" -R "{N.Reference_Fasta}"'
	BashSubprocess(Name = f'HaplotypeCalling.Contig.{Contig}.HaplotypeCaller', Command = CommandHaplotypeCaller)
	BashSubprocess(Name = f'HaplotypeCalling.Contig.{Contig}.GenotypeGVCFs', Command = CommandGenotypeGVCFs)
	return { 'VCF': OutputVCF, 'gVCF': OutputGVCF }

def HaplotypeCalling(**kwargs):
	logging.info('* entry point *')
	N = RenderParameters(kwargs)
	with tempfile.TemporaryDirectory() as TempDir:
		with Threading('ContigHaplotypeCalling', N.Threads) as pool:
			Shards = pool.map(functools.partial(ContigHaplotypeCalling, Input_BAM = N.Input_BAM, Temp_Directory = TempDir, Reference_Fasta = N.Reference_Fasta), N.Active_Contigs)
			TaggedShardsVCF = MultipleTags(Tag = '-I', List = [item['VCF'] for item in Shards])
			TaggedShardsGVCF = MultipleTags(Tag = '-I', List = [item['gVCF'] for item in Shards])
		CommandMergeVcfs = f'"{GATK_PATH}" --java-options "{CONFIG_JAVA_OPTIONS["MergeVcfs"]}" MergeVcfs --CREATE_INDEX true --VERBOSITY ERROR --CREATE_MD5_FILE true {TaggedShardsVCF} -O "{N.Output_VCF}"'
		CommandMergeGVcfs = f'"{GATK_PATH}" --java-options "{CONFIG_JAVA_OPTIONS["MergeVcfs"]}" MergeVcfs --CREATE_INDEX true --VERBOSITY ERROR --CREATE_MD5_FILE true {TaggedShardsGVCF} -O "{N.Output_gVCF}"'
		BashSubprocess(Name = 'HaplotypeCalling.MergeVCF', Command = CommandMergeVcfs)
		BashSubprocess(Name = 'HaplotypeCalling.MergeGVCF', Command = CommandMergeGVcfs)
