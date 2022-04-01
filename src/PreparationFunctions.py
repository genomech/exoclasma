from .SharedFunctions import *


## ------======| REFSEQ PREP |======------

def CreateGenomeInfo(GenomeName, Dir):
	logging.info(RenderParameters('*', '*'))
	ConfigJson = {
		'NAME':       '',
		'FASTA':      f'{GenomeName}.fasta',
		'CHROMSIZES': f'{GenomeName}.chromsizes.txt',
		'FAI':        f'{GenomeName}.fasta.fai',
		'BED':        f'{GenomeName}.bed',
		'DICT':       f'{GenomeName}.dict'
	}
	ConfigJson = { Key: os.path.join(Dir, Value) for Key, Value in ConfigJson.items() }
	ConfigJson['RS'] = { Name: f'{GenomeName}.restriction_sites.{Name}.txt' for Name in CONFIG_RESTRICTION_ENZYMES['Enzymes'].keys() }
	ConfigJson['RS'] = { Name: os.path.join(Dir, FN) for Name, FN in ConfigJson['RS'].items() }
	ConfigJson['NAME'] = GenomeName
	return ConfigJson

def RefseqPreparation(**kwargs):
	logging.info(RenderParameters('*', '*'))
	for Key, Value in kwargs.items(): logging.info(RenderParameters(Key, Value))
	N = types.SimpleNamespace(**kwargs)
	OutputDir = os.path.join(N.Parent_Dir, N.Genome_Name)
	os.mkdir(OutputDir)
	logging.info(RenderParameters('Directory created', OutputDir))
	Fasta = SeqIO.parse(OpenAnyway(N.Fasta_Path, 'rt'), 'fasta')
	logging.info(RenderParameters('FASTA opened', N.Fasta_Path))
	SearchQueries = { name: re.compile(sequences) for name, sequences in CONFIG_RESTRICTION_ENZYMES['Enzymes'].items() }
	with tempfile.TemporaryDirectory() as TempDir:
		FNdict = CreateGenomeInfo(N.Genome_Name, TempDir)
		with open(FNdict['FASTA'], 'w') as NewFasta, open(FNdict['CHROMSIZES'], 'w') as ChromSizes:
			for contig in Fasta:
				contig.name = re.sub('[^\w\.]', '_', contig.name)
				Seq = contig.seq.__str__()
				SeqLength = len(Seq)
				SeqIO.write([contig], NewFasta, 'fasta')
				ChromSizes.write(f'{contig.name}\t{SeqLength}\n')
				for enzyme, query in SearchQueries.items():
					FileWrapper = open(FNdict['RS'][enzyme], 'a')
					FileWrapper.write(' '.join([contig.name] + [str(match.start() + 1) for match in query.finditer(Seq)] + [str(SeqLength)]) + '\n')
					FileWrapper.close()
				logging.info(RenderParameters('Contig ready', contig.name))
		for Key, Value in {'File ready': 'Fasta', 'File ready': 'ChromSizes', 'Files ready': 'RestrictionSites'}.items(): logging.info(RenderParameters(Key, Value))
		CommandSamtoolsIndex = f'samtools faidx "{FNdict["FASTA"]}"'
		CommandBwaIndex = f'bwa index "{FNdict["FASTA"]}"'
		CommandGATKIndex = f'"{GATK_PATH}" CreateSequenceDictionary --VERBOSITY ERROR -R "{FNdict["FASTA"]}"'
		CommandGenomeBED = f'awk \'BEGIN {{ FS = "\\t" }}; {{ print $1 FS "0" FS $2 }}\' "{FNdict["FAI"]}" > "{FNdict["BED"]}"'
		CommandCopy = f'mv {os.path.join(TempDir, "*")} "{OutputDir}"'
		BashSubprocess(Name = 'RefseqPreparation.SAMtoolsIndex', Command = CommandSamtoolsIndex)
		BashSubprocess(Name = 'RefseqPreparation.BWAIndex', Command = CommandBwaIndex)
		BashSubprocess(Name = 'RefseqPreparation.GATKIndex', Command = CommandGATKIndex)
		BashSubprocess(Name = 'RefseqPreparation.GenomeBED', Command = CommandGenomeBED)
		JsonFN = os.path.join(TempDir, f'{N.Genome_Name}.info.json')
		ConfigJson = CreateGenomeInfo(N.GenomeName, OutputDir)
		SaveJSON(ConfigJson, JsonFN)
		BashSubprocess(Name = 'RefseqPreparation.Copy', Command = CommandCopy)


## ------======| CAPTURE PREP |======------

def CreateCaptureInfo(CaptureName, Dir):
	logging.info(RenderParameters('*', '*'))
	ConfigJson = {
		'NAME':     '',
		'CAP':      f'{CaptureName}.capture.bed',
		'NOTCAP':   f'{CaptureName}.not_capture.bed'
	}
	ConfigJson = { Key: os.path.join(Dir, Value) for Key, Value in ConfigJson.items() }
	ConfigJson['NAME'] = CaptureName
	return ConfigJson

def CapturePreparation(**kwargs):
	logging.info(RenderParameters('*', '*'))
	for Key, Value in kwargs.items(): logging.info(RenderParameters(Key, Value))
	N = types.SimpleNamespace(**kwargs)
	OutputDir = os.path.join(N.Parent_Dir, N.Capture_Name)
	os.mkdir(OutputDir)
	logging.info(RenderParameters('Directory created', OutputDir))
	GenomeInfo = LoadJSON(N.Genome_Info_JSON)
	with tempfile.TemporaryDirectory() as TempDir:
		FNdict = CreateCaptureInfo(N.Capture_Name, TempDir)
		CommandFilterAndSort = f'set -o pipefail; bedtools intersect -a "{GenomeInfo["BED"]}" -b "{N.Input_BED}" | bedtools sort -faidx "{GenomeInfo["FAI"]}" | sed -e \'s/$/\\t\\./\' > "{FNdict["CAP"]}"'
		CommandNotCapture = f'bedtools subtract -a "{GenomeInfo["BED"]}" -b "{FNdict["CAP"]}" | sed -e \'s/$/\\t\\./\' > "{FNdict["NOTCAP"]}"'
		CommandCopy = f'mv {os.path.join(TempDir, "*")} "{OutputDir}"'
		BashSubprocess(Name = 'CapturePreparation.FilterSort', Command = CommandFilterAndSort)
		BashSubprocess(Name = 'CapturePreparation.NotCapture', Command = CommandNotCapture)
		JsonFN = os.path.join(TempDir, f'{N.Capture_Name}.info.json')
		ConfigJson = CreateCaptureInfo(N.Capture_Name, OutputDir)
		SaveJSON(ConfigJson, JsonFN)
		BashSubprocess(Name = 'CapturePreparation.Copy', Command = CommandCopy)
