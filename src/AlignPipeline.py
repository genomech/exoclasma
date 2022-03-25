from .SharedFunctions import *
from .AlignFunctions import *
from .Statistics import *


# ------======| PROTOCOL PREPARATION & BACKUP |======------

def GenerateAlignFileNames(Unit, Meta):
	ID = Unit['ID']
	OutputDir = os.path.join(Meta['PoolDir'], ID)
	FileNames = {
		"Log":             f"{ID}.log",
		"RecalBAM":        f"{ID}.bam",
		"VCF":             f"{ID}.vcf.gz",
		"gVCF":            f"{ID}.g.vcf.gz",
		"FullStats":       f"{ID}.stats.json",
		"PrimaryBAM":      f"_temp.{ID}.primary.bam",
		"FlagStat":        f"_temp.{ID}.flagstat.json",
		"DuplessBAM":      f"_temp.{ID}.dupless.bam",
		"DuplessQSBAM":    f"_temp.{ID}.duplessqs.bam",
		"DuplessMetrics":  f"_temp.{ID}.md_metrics.txt",
		"CoverageStats":   f"_temp.{ID}.coverage.json",
		"Contigs":         f"_temp.{ID}.contigs.json"
	}
	FileNames = { Key: os.path.join(OutputDir, Value) for Key, Value in FileNames.items() }
	FileNames["Cutadapt"] = { str(index): {
		'FastQC':          f"{ID}.{index}.trimmed_fastqc.html",
		'R1':              f"_temp.{ID}.{index}.R1.fastq.gz",
		'R2':              f"_temp.{ID}.{index}.R2.fastq.gz",
		'Stats':           f"_temp.{ID}.{index}.coverage.tsv",
		'ActiveContigs':   f"_temp.{ID}.{index}.contigs.txt"
		} for index in range(len(Unit['Input'])) if Unit['Input'][index]['Type'] == 'fastq' }
	FileNames["Cutadapt"] = { Key: { Key2: os.path.join(OutputDir, Value2) for Key2, Value2 in Value.items() } for Key, Value in FileNames["Cutadapt"].items() }
	FileNames["OutputDir"] = OutputDir
	return FileNames

def AlignProtocolAndBackup(UnitsFile):
	Backup = f"{UnitsFile}.align.backup"
	if os.path.exists(Backup) and os.path.isfile(Backup):
		Protocol = LoadJSON(Backup)
		logging.warning(f'Resume process from backup "{Backup}"')
	else:
		Protocol = LoadJSON(UnitsFile)
		Protocol['Meta']['GenomeInfo'] = LoadJSON(Protocol['Meta']['GenomeInfo'])
		Protocol['Meta']['CaptureInfo'] = LoadJSON(Protocol['Meta']['CaptureInfo'])
		for index in range(len(Protocol['Units'])): 
			Protocol['Units'][index]['FileNames'] = GenerateAlignFileNames(Protocol['Units'][index], Protocol['Meta'])
		Protocol['Backup'] = { 'Possible': True, 'FN': Backup }
	try:
		SaveJSON(Protocol, Backup)
	except OSError:
		logging.warning(f'Backup file "{Backup}" cannot be created. Check current user permissions.')
		Protocol['Backup']['Possible'] = False
	return Protocol

# ------======| STAGES |======------

def MakeDirStage(Unit, Meta): os.mkdir(Unit['FileNames']["OutputDir"])

def CutadaptStage(Unit, Meta):
	for index, item in enumerate(Unit["Input"]):
		index = str(index)
		if (item["Type"] == "fastq") and (item["Adapter"] is not None):
			Cutadapt(
				InputR1 = item['R1'],
				InputR2 = item['R2'],
				OutputR1 = Unit['FileNames']['Cutadapt'][index]['R1'],
				OutputR2 = Unit['FileNames']['Cutadapt'][index]['R2'],
				Adapter = CONFIG_ADAPTERS[item["Adapter"]],
				ReportTXT = Unit['FileNames']['Cutadapt'][index]['Stats'],
				Threads = Meta["Threads"]
			)
			FastQC(
				InputFastQ = Unit['FileNames']['Cutadapt'][index]['R1'],
				OutputHTML = Unit['FileNames']['Cutadapt'][index]['FastQC'],
				Size = Meta["FastQSampleSize"],
				Threads = Meta["Threads"]
			)

def BWAStage(Unit, Meta):
	with tempfile.TemporaryDirectory() as TempDir:
		Shards = []
		for index, item in enumerate(Unit["Input"]):
			index = str(index)
			OutputBAM = os.path.join(TempDir, f"temp_{index}.bam")
			if item["Type"] == "fastq":
				RGtag = RGTag(item['Instrument'], item['Lane'], item['Platform'], item['Barcode'], item['Sample'], item['Library'])
				if item["Adapter"] is not None: 
					R1 = Unit['FileNames']['Cutadapt'][index]['R1']
					R2 = Unit['FileNames']['Cutadapt'][index]['R2']
				else: 
					R1 = item['R1']
					R2 = item['R2']
				BWA(
					InputR1 = R1,
					InputR2 = R2,
					Reference = Meta["GenomeInfo"]['FASTA'],
					RGHeader = RGtag,
					OutputBAM = OutputBAM,
					ActiveContigsTXT = Unit['FileNames']['Cutadapt'][index]['ActiveContigs'],
					Threads = Meta["Threads"]
				)
				Shards += [ OutputBAM ]
			if item["Type"] == "bam": pass
		if len(Shards) > 1:
			MergeBAMs(
				BAMs = Shards,
				OutputBAM = Unit['FileNames']["PrimaryBAM"],
				SortOrder = "queryname"
			)
		else: 
			BashSubprocess(
				Name = f'BWAStage.CopyBAM',
				Command = f'cp "{Shards[0]}" "{Unit["FileNames"]["PrimaryBAM"]}"'
			)
		FlagStat(InputBAM = Unit['FileNames']["PrimaryBAM"], OutputJSON = Unit['FileNames']["FlagStat"], Threads = Meta["Threads"])

def GetActiveContigs(Unit, Meta):
	Data = [ LoadContigList(item['ActiveContigs']) for item in Unit['FileNames']['Cutadapt'].values() ]
	Data = list(set([item for sublist in Data for item in sublist]))
	Chroms = pandas.read_csv(Meta['GenomeInfo']['CHROMSIZES'], sep='\t', header=None)[0].to_list()
	SortedContigs = [item for item in Chroms if item in Data]
	SaveJSON(SortedContigs, Unit["FileNames"]["Contigs"])
	ContigsString = ', '.join(SortedContigs)
	logging.info(f'Active contigs: {ContigsString}')

def MarkDuplicatesStage(Unit, Meta):
	MarkDuplicates(
		InputBAM = Unit['FileNames']["PrimaryBAM"],
		OutputBAM = Unit['FileNames']["DuplessBAM"],
		QuerySortedBAM = Unit['FileNames']["DuplessQSBAM"]
		MetricsTXT = Unit['FileNames']["DuplessMetrics"]
	)

def CoverageStatsStage(Unit, Meta):
	CoverageStats(
		Name = Unit['ID'],
		FinalBAM = Unit['FileNames']["DuplessBAM"],
		StatsTXT = Unit['FileNames']['CoverageStats'],
		CaptureInfo = Meta['CaptureInfo'],
		GenomeInfo = Meta['GenomeInfo']
	)

def BaseRecalibrationStage(Unit, Meta):
	BaseRecalibration(
		InputBAM = Unit['FileNames']["DuplessBAM"],
		OutputBAM = Unit['FileNames']["RecalBAM"],
		dbSNP = None,
		Reference = Meta["GenomeInfo"]['FASTA'],
		ActiveContigs = Unit['ActiveContigs'],
		Threads = Meta["Threads"]
	) 

def DuplessAsFinal(Unit):
	BashSubprocess(
		Name = f'DuplessAsFinal.Move',
		Command = f'mv "{Unit["FileNames"]["DuplessBAM"]}" "{Unit["FileNames"]["RecalBAM"]}"'
			)

def HaplotypeCallingStage(Unit, Meta):
	HaplotypeCalling(
		InputBAM = Unit['FileNames']["RecalBAM"],
		OutputVCF = Unit['FileNames']["VCF"],
		OutputGVCF = Unit['FileNames']["gVCF"],
		Reference = Meta["GenomeInfo"]['FASTA'],
		ActiveContigs = Unit['ActiveContigs'],
		Threads = Meta["Threads"]
	)

def StatsSummaryStage(Unit, Meta):
	Result = {
		'RefseqName':      Meta['GenomeInfo']['NAME'],
		'CaptureName':     Meta['CaptureInfo']['NAME'],
		'Cutadapt':        { index: CutadaptStat(item['Stats']) for index, item in Unit['FileNames']['Cutadapt'].items() },
		'FlagStat':        LoadFlagStat(Unit['FileNames']['FlagStat']),
		'MarkDuplicates':  LoadMarkDuplicatesStat(Unit['FileNames']['DuplessMetrics']),
		'CoverageStats':   LoadCoverageStats(Unit['FileNames']['CoverageStats'])
		}
	SaveJSON(Result, Unit['FileNames']['FullStats'])

# ------======| ALIGN PIPELINE |======------

def AlignPipeline(UnitsFile, Verbosity = logging.INFO):
	Protocol = AlignProtocolAndBackup(UnitsFile)
	for UnitIndex in range(len(Protocol["Units"])):
		Parameters = { 'Unit': Protocol["Units"][UnitIndex], 'Meta': Protocol['Meta'] }
		def StageAndBackup(Stage, Func, Parameters):
			nonlocal Protocol
			nonlocal UnitIndex
			if Stage not in Protocol["Units"][UnitIndex]["Stage"]:
				Func(**Parameters)
				Protocol["Units"][UnitIndex]["Stage"].append(Stage)
				if Protocol['Backup']['Possible']: SaveJSON(Protocol, Protocol['Backup']['FN'])
		StartTime = time.time()
		StageAndBackup(Stage = "MakeDir", Func = MakeDirStage, Parameters = Parameters)
		ConfigureLogger(Protocol["Units"][UnitIndex]['FileNames']['Log'], Verbosity)
		StageAndBackup(Stage = "Cutadapt", Func = CutadaptStage, Parameters = Parameters)
		StageAndBackup(Stage = "BWA", Func = BWAStage, Parameters = Parameters)
		StageAndBackup(Stage = "GetActiveContigs", Func = GetActiveContigs, Parameters = Parameters)
		Protocol["Units"][UnitIndex]['ActiveContigs'] = LoadJSON(Protocol["Units"][UnitIndex]["FileNames"]["Contigs"])
		Parameters = { 'Unit': Protocol["Units"][UnitIndex], 'Meta': Protocol['Meta'] }
		StageAndBackup(Stage = "MarkDuplicates", Func = MarkDuplicatesStage, Parameters = Parameters)
		StageAndBackup(Stage = "CoverageStats", Func = CoverageStatsStage, Parameters = Parameters)
		if Protocol['Meta']['Call']:
			StageAndBackup(Stage = "BaseRecalibration", Func = BaseRecalibrationStage, Parameters = Parameters)
			StageAndBackup(Stage = "HaplotypeCalling", Func = HaplotypeCallingStage, Parameters = Parameters)
		else: DuplessAsFinal(Parameters['Unit'])
		if Protocol['Meta']['HiC']:
			# TODO
		StageAndBackup(Stage = "StatsSummaryStage", Func = StatsSummaryStage, Parameters = Parameters)
		if Protocol['Meta']['RemoveTempFiles']:
			for t in glob.glob(os.path.join(Protocol["Units"][UnitIndex]['FileNames']["OutputDir"], '_temp.*')): os.remove(t)
			logging.info(f'Temp files removed')
		logging.info(f'Unit: "{Protocol["Units"][UnitIndex]["ID"]}"; Summary time: {Timestamp(StartTime)}')
	
	os.remove(Protocol['Backup']['FN'])
