from .SharedFunctions import *
from .AlignFunctions import *
from .Statistics import *
from .Jooser import *


# ------======| PROTOCOL PREPARATION & BACKUP |======------

def GenerateAlignFileNames(Unit):
	logging.info(RenderParameters('*', '*'))
	ID = Unit['ID']
	OutputDir = os.path.join(Unit['PoolDir'], ID)
	FileNames = {
		'Log':              f'{ID}.log',
		'RecalBAM':         f'{ID}.bam',
		'VCF':              f'{ID}.vcf.gz',
		'gVCF':             f'{ID}.g.vcf.gz',
		'MergedNoDups':     f'{ID}.merged_nodups.txt.gz',
		'Inter30':          f'{ID}.inter30.hic',
		'FullStats':        f'{ID}.stats.json',
		'PrimaryBAM':       f'_temp.{ID}.primary.bam',
		'FlagStat':         f'_temp.{ID}.flagstat.json',
		'DuplessBAM':       f'_temp.{ID}.dupless.bam',
		'DuplessQSBAM':     f'_temp.{ID}.duplessqs.bam',
		'DuplessMetrics':   f'_temp.{ID}.md_metrics.txt',
		'CoverageStats':    f'_temp.{ID}.coverage.json',
		'Contigs':          f'_temp.{ID}.contigs.json',
		'JooserStats':      f'_temp.{ID}.jstats.json'
	}
	FileNames = { Key: os.path.join(OutputDir, Value) for Key, Value in FileNames.items() }
	FileNames['Cutadapt'] = { str(index): {
		'R1':               f'{ID}.{index}.R1.fastq.gz',
		'R2':               f'{ID}.{index}.R2.fastq.gz',
		'Unpaired':         f'{ID}.{index}.unpaired.fastq.gz',
		'FastQC':           f'{ID}.{index}.fastqc.html',
		'CutadaptR1':       f'_temp.{ID}.{index}.R1.cutadapt.fastq.gz',
		'CutadaptR2':       f'_temp.{ID}.{index}.R2.cutadapt.fastq.gz',
		'CutadaptUnpaired': f'_temp.{ID}.{index}.unpaired.cutadapt.fastq.gz',
		'Stats':            f'_temp.{ID}.{index}.cutstats.tsv'
		} for index in range(len(Unit['Input'])) }
	FileNames['Cutadapt'] = { Key: { Key2: os.path.join(OutputDir, Value2) for Key2, Value2 in Value.items() } for Key, Value in FileNames['Cutadapt'].items() }
	FileNames['OutputDir'] = OutputDir
	return FileNames

def AlignProtocolAndBackup(UnitsFile):
	logging.info(RenderParameters('*', '*'))
	global GLOBAL_BACKUP
	Backup = os.path.join(os.path.dirname(os.path.realpath(UnitsFile)), f'{UnitsFile}.align.backup')
	GLOBAL_BACKUP = { 'Possible': True, 'FN': Backup }
	if os.path.exists(Backup) and os.path.isfile(Backup):
		Protocol = LoadJSON(Backup)
		logging.warning(RenderParameters('Resume process', Backup))
	else:
		Protocol = LoadJSON(UnitsFile)
		for index in range(len(Protocol)): 
			Protocol[index]['Reference']['GenomeInfo'] = LoadJSON(Protocol[index]['Reference']['GenomeInfo'])
			Protocol[index]['Reference']['CaptureInfo'] = LoadJSON(Protocol[index]['Reference']['CaptureInfo'])
			Protocol[index]['Output'] = GenerateAlignFileNames(Protocol[index])
	try:
		SaveJSON(Protocol, Backup)
	except OSError:
		logging.warning(f'Backup file "{Backup}" cannot be saved. Check current user permissions.')
		GLOBAL_BACKUP['Possible'] = False
	return Protocol


# ------======| STAGES |======------

def MakeDirStage(Unit):
	logging.info(RenderParameters('*', '*'))
	os.mkdir(Unit['Output']['OutputDir'])

def SraToFastQStage(Unit):
	logging.info(RenderParameters('*', '*'))
	Info = {}
	for index, item in enumerate(Unit['Input']):
		index = str(index)
		if item['Type'] == 'sra':
			Info[int(index)] = SraToFastq(
				SRA_Dataset = item['Files']['SRA'],
				SRA_FastQ_R1 = Unit['Output']['Cutadapt'][index]['R1'],
				SRA_FastQ_R2 = Unit['Output']['Cutadapt'][index]['R2'],
				SRA_FastQ_Unpaired = Unit['Output']['Cutadapt'][index]['Unpaired']
				)
	return Info

def CutadaptStage(Unit):
	logging.info(RenderParameters('*', '*'))
	for index, item in enumerate(Unit['Input']):
		index = str(index)
		if item['Adapter'] is not None:
			if item['Files']['Unpaired'] is not None:
				Cutadapt(
					Mode = 'Single-end',
					Input_FastQ = item['Files']['Unpaired'],
					Output_FastQ = Unit['Output']['Cutadapt'][index]['CutadaptUnpaired'],
					Adapter = item['Adapter'],
					Cutadapt_Report = Unit['Output']['Cutadapt'][index]['Stats'],
					Threads = Unit['Config']['Threads']
					)
				FileToAnalyze = Unit['Output']['Cutadapt'][index]['CutadaptUnpaired']
			else: Unit['Output']['Cutadapt'][index]['CutadaptUnpaired'] = None
			if (item['Files']['R1'] is not None) and (item['Files']['R2'] is not None):
				Cutadapt(
					Mode = 'Paired-end',
					Input_FastQ_R1 = item['Files']['R1'],
					Input_FastQ_R2 = item['Files']['R2'],
					Output_FastQ_R1 = Unit['Output']['Cutadapt'][index]['CutadaptR1'],
					Output_FastQ_R2 = Unit['Output']['Cutadapt'][index]['CutadaptR2'],
					Adapter = item['Adapter'],
					Cutadapt_Report = Unit['Output']['Cutadapt'][index]['Stats'],
					Threads = Unit['Config']['Threads']
					)
				FileToAnalyze = Unit['Output']['Cutadapt'][index]['CutadaptR1']
			else:
				Unit['Output']['Cutadapt'][index]['CutadaptR1'] = None
				Unit['Output']['Cutadapt'][index]['CutadaptR2'] = None
			FastQC(
				Input_FastQ = FileToAnalyze,
				Output_HTML = Unit['Output']['Cutadapt'][index]['FastQC'],
				Subsample_Size = Unit['Config']['FastQSampleSize'],
				Threads = Unit['Config']['Threads']
				)

def BWAStage(Unit):
	logging.info(RenderParameters('*', '*'))
	with tempfile.TemporaryDirectory() as TempDir:
		Shards = []
		for index, item in enumerate(Unit['Input']):
			index = str(index)
			RGtag = RGTag(**(item['RG']))
			if item['Adapter'] is not None: 
				R1 = Unit['Output']['Cutadapt'][index]['CutadaptR1']
				R2 = Unit['Output']['Cutadapt'][index]['CutadaptR2']
				Unpaired = Unit['Output']['Cutadapt'][index]['CutadaptUnpaired']
			else: 
				R1 = item['Files']['R1']
				R2 = item['Files']['R2']
				Unpaired = item['Files']['Unpaired']
			if Unpaired is not None:
				OutputBAM = os.path.join(TempDir, f'temp_{index}_unpaired.bam')
				BWA(
					Mode = 'Single-end',
					Input_FastQ = Unpaired,
					Reference_Fasta = Unit['Reference']['GenomeInfo']['FASTA'],
					RG_Header = RGtag,
					Output_BAM = OutputBAM,
					Threads = Unit['Config']['Threads']
					)
				Shards.append(OutputBAM)
			if (R1 is not None) and (R2 is not None):
				OutputBAM = os.path.join(TempDir, f'temp_{index}_paired.bam')
				BWA(
					Mode = 'Paired-end',
					Input_FastQ_R1 = R1,
					Input_FastQ_R2 = R2,
					Reference_Fasta = Unit['Reference']['GenomeInfo']['FASTA'],
					RG_Header = RGtag,
					Output_BAM = OutputBAM,
					Threads = Unit['Config']['Threads']
					)
				Shards.append(OutputBAM)
		if len(Shards) > 1:
			MergeSamFiles(
				BAM_List = Shards,
				Output_BAM = Unit['Output']['PrimaryBAM'],
				Sort_Order = 'queryname'
				)
		else:
			CommandMove = f'mv "{Shards[0]}" "{Unit["Output"]["PrimaryBAM"]}"'
			BashSubprocess(Name = f'BWAStage.MoveBAM', Command = CommandMove)
		FlagStat(Input_BAM = Unit['Output']['PrimaryBAM'], Samtools_Flagstats = Unit['Output']['FlagStat'], Threads = Unit['Config']['Threads'])

def GetActiveContigs(Unit):
	logging.info(RenderParameters('*', '*'))
	Data = [item for item in ''.join(open(Unit['Output']['Contigs'], 'rt').readlines())[:-1].split('\n') if item != '*']
	Chroms = pandas.read_csv(Unit['Reference']['GenomeInfo']['CHROMSIZES'], sep = '\t', header = None)[0].to_list()
	SortedContigs = [item for item in Chroms if item in Data]
	SaveJSON(SortedContigs, Unit['Output']['Contigs'])
	ContigsString = ', '.join(SortedContigs)
	logging.info(RenderParameters('Active Contigs', ContigsString))

def MarkDuplicatesStage(Unit):
	logging.info(RenderParameters('*', '*'))
	MarkDuplicates(
		Input_BAM = Unit['Output']['PrimaryBAM'],
		Output_BAM = Unit['Output']['DuplessBAM'],
		Query_Sorted_BAM = Unit['Output']['DuplessQSBAM'],
		MarkDup_Metrics = Unit['Output']['DuplessMetrics'],
		Active_Contigs_File = Unit['Output']['Contigs']
		)

def CoverageStatsStage(Unit):
	logging.info(RenderParameters('*', '*'))
	CoverageStats(
		Input_BAM = Unit['Output']['DuplessBAM'],
		Coverage_Stats = Unit['Output']['CoverageStats'],
		Capture_BED = Unit['Reference']['CaptureInfo']['CAP'],
		NotCapture_BED = Unit['Reference']['CaptureInfo']['NOTCAP'],
		Reference_FAI = Unit['Reference']['GenomeInfo']['FAI']
		)

def BaseRecalibrationStage(Unit):
	logging.info(RenderParameters('*', '*'))
	BaseRecalibration(
		Input_BAM = Unit['Output']['DuplessBAM'],
		Output_BAM = Unit['Output']['RecalBAM'],
		dbSNP_Known_Sites = '/Data/Tools/exoclasma/databases/db/dbsnp151_BaseRecalibrator.vcf.gz',
		Reference_Fasta = Unit['Reference']['GenomeInfo']['FASTA'],
		Active_Contigs = Unit['ActiveContigs'],
		Threads = Unit['Config']['Threads']
		) 

def DuplessAsFinal(Unit):
	logging.info(RenderParameters('*', '*'))
	Command = f'mv "{Unit["FileNames"]["DuplessBAM"]}" "{Unit["FileNames"]["RecalBAM"]}"'
	BashSubprocess(Name = 'DuplessAsFinal', Command = Command)

def HaplotypeCallingStage(Unit):
	logging.info(RenderParameters('*', '*'))
	HaplotypeCalling(
		Input_BAM = Unit['Output']['RecalBAM'],
		Output_VCF = Unit['Output']['VCF'],
		Output_gVCF = Unit['Output']['gVCF'],
		Reference_Fasta = Unit['Reference']['GenomeInfo']['FASTA'],
		Active_Contigs = Unit['ActiveContigs'],
		Threads = Unit['Config']['Threads']
		)

def MergedNoDupsStage(Unit):
	logging.info(RenderParameters('*', '*'))
	JooserFunc(
		Input_BAM = Unit['Output']['DuplessQSBAM'],
		MergedNoDups_File = Unit['Output']['MergedNoDups'],
		Jooser_Stats = Unit['Output']['JooserStats'],
		Restriction_Site_Map = None if Unit['RS'] is None else Unit['Reference']['GenomeInfo']['RS'][Unit['RS']],
		Min_MAPQ = Unit['Config']['MinMAPQ']
		)

def JuicerToolsStage(Unit):
	logging.info(RenderParameters('*', '*'))
	JuicerTools(
		MergedNoDups_File = Unit['Output']['MergedNoDups'],
		Output_HIC_File = Unit['Output']['Inter30'],
		ChromSizes_File = Unit['Reference']['GenomeInfo']['CHROMSIZES'],
		Restriction_Site_Map = None if Unit['RS'] is None else Unit['Reference']['GenomeInfo']['RS'][Unit['RS']],
		Threads = Unit['Config']['Threads']
		)

def StatsSummaryStage(Unit):
	logging.info(RenderParameters('*', '*'))
	Result = {
		'RefseqName':      Unit['Reference']['GenomeInfo']['NAME'],
		'CaptureName':     Unit['Reference']['CaptureInfo']['NAME'],
		'FlagStat':        LoadJSON(Unit['Output']['FlagStat']),
		'MarkDuplicates':  LoadMarkDuplicatesStat(Unit['Output']['DuplessMetrics']),
		'CoverageStats':   LoadJSON(Unit['Output']['CoverageStats'])
		}
	if any([item['Adapter'] is not None for item in Unit['Input']]):
		Result['Cutadapt'] = { index: CutadaptStat(item['Stats']) for index, item in Unit['Output']['Cutadapt'].items() if os.path.exists(item['Stats']) }
	if Unit['Config']['HiC']: Result['JooserStats'] = LoadJSON(Unit['Output']['JooserStats'])
	SaveJSON(Result, Unit['Output']['FullStats'])

# ------======| ALIGN PIPELINE |======------

def AlignPipeline(UnitsFile, Verbosity = logging.INFO):
	Protocol = AlignProtocolAndBackup(UnitsFile)
	for UnitIndex in range(len(Protocol)):
		def StageAndBackup(Stage, Func):
			nonlocal Protocol
			nonlocal UnitIndex
			global GLOBAL_BACKUP
			if Stage not in Protocol[UnitIndex]['Stage']:
				TechInfo = Func(Protocol[UnitIndex])
				Protocol[UnitIndex]['Stage'].append(Stage)
				if GLOBAL_BACKUP['Possible']: SaveJSON(Protocol, GLOBAL_BACKUP['FN'])
				return TechInfo
			return {}
		StartTime = time.time()
		StageAndBackup(Stage = 'MakeDir', Func = MakeDirStage)
		ConfigureLogger(Protocol[UnitIndex]['Output']['Log'], Verbosity)
		SRAInfo = StageAndBackup(Stage = 'SraToFastQ', Func = SraToFastQStage)
		for Index, Item in SRAInfo.items(): Protocol[UnitIndex]['Input'][Index]['Files'] = Item
		if GLOBAL_BACKUP['Possible']: SaveJSON(Protocol, GLOBAL_BACKUP['FN'])
		StageAndBackup(Stage = 'Cutadapt', Func = CutadaptStage)
		StageAndBackup(Stage = 'BWA', Func = BWAStage)
		StageAndBackup(Stage = 'MarkDuplicates', Func = MarkDuplicatesStage)
		StageAndBackup(Stage = 'GetActiveContigs', Func = GetActiveContigs)
		Protocol[UnitIndex]['ActiveContigs'] = LoadJSON(Protocol[UnitIndex]['Output']['Contigs'])
		StageAndBackup(Stage = 'CoverageStats', Func = CoverageStatsStage)
		if Protocol[UnitIndex]['Config']['Call']:
			StageAndBackup(Stage = 'BaseRecalibration', Func = BaseRecalibrationStage)
			StageAndBackup(Stage = 'HaplotypeCalling', Func = HaplotypeCallingStage)
		else: StageAndBackup(Stage = 'DuplessAsFinal', Func = DuplessAsFinal)
		if Protocol[UnitIndex]['Config']['HiC']:
			StageAndBackup(Stage = 'MergedNoDups', Func = MergedNoDupsStage)
			StageAndBackup(Stage = 'JuicerTools', Func = JuicerToolsStage)
		StageAndBackup(Stage = 'StatsSummaryStage', Func = StatsSummaryStage)
		if Protocol[UnitIndex]['Config']['RemoveTempFiles']:
			for t in glob.glob(os.path.join(Protocol[UnitIndex]['Output']['OutputDir'], '_temp.*')): os.remove(t)
			logging.info(f'Temp files removed')
		logging.info(f'Summary time: {Timestamp(StartTime)}')
	os.remove(GLOBAL_BACKUP['FN'])
