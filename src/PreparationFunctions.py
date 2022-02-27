from src.SharedFunctions import *

## ------======| REFSEQ PREP |======------

def SeriesReplace(String, Replacements):
	Rep = dict((re.escape(key), value) for key, value in Replacements.items()) 
	return re.compile("|".join(Rep.keys())).sub(lambda x: Rep[re.escape(x.group(0))], String)

def CreateGenomeInfo(GenomeName, Dir):
	ConfigJson = {
		'NAME':       GenomeName,
		'FASTA':      f'{GenomeName}.fasta',
		'CHROMSIZES': f'{GenomeName}.chromsizes.txt',
		'FAI':        f'{GenomeName}.fasta.fai',
		'BED':        f'{GenomeName}.bed',
		'DICT':       f'{GenomeName}.dict'
	}
	ConfigJson = { Key: os.path.join(Dir, Value) for Key, Value in ConfigJson }
	ConfigJson['RS'] = { Name: f'{GenomeName}.restriction_sites.{Name}.txt' for Name in CONFIG_RESTRICTION_ENZYMES["Enzymes"].keys() }
	ConfigJson['RS'] = { Name: os.path.join(Dir, FN) for Name, FN in ConfigJson['RS'].items() }
	return ConfigJson

def RefseqPreparation(FastaPath, GenomeName, ParentDir, Threads = multiprocessing.cpu_count()):
	OutputDir = os.path.join(ParentDir, GenomeName)
	os.mkdir(OutputDir)
	logging.info(f'Directory created: "{OutputDir}"')
	Fasta = SeqIO.parse(OpenAnyway(FastaPath, "rt"), "fasta")
	logging.info(f'FASTA opened: "{FastaPath}"')
	SearchQueries = { name: re.compile("|".join([f"({SeriesReplace(item, CONFIG_RESTRICTION_ENZYMES['WildCards'])})" for item in sequences])) for name, sequences in CONFIG_RESTRICTION_ENZYMES["Enzymes"].items() }
	with tempfile.TemporaryDirectory() as TempDir:
		FNdict = CreateGenomeInfo(GenomeName, TempDir)
		with open(FNdict['FASTA'], "w") as NewFasta, open(FNdict['CHROMSIZES'], "w") as ChromSizes:
			for contig in Fasta:
				contig.name = re.sub("[^\w\.]", "_", contig.name)
				Seq = contig.seq.__str__()
				SeqLength = len(Seq)
				SeqIO.write([contig], NewFasta, 'fasta')
				ChromSizes.write(f"{contig.name}\t{SeqLength}\n")
				for enzyme, query in SearchQueries.items():
					FileWrapper = open(FNdict['RS'][enzyme], "a")
					FileWrapper.write(" ".join([contig.name] + [str(match.start() + 1) for match in query.finditer(Seq)] + [str(SeqLength)]) + "\n")
					FileWrapper.close()
				logging.info(f'Contig "{contig.name}" ready')
		logging.info(f"New FASTA, ChromSizes and RestrictionSites are ready")
		BashSubprocess(Name = f"SAMtools Index", Command = f'samtools faidx -@ {Threads} "{FNdict["FASTA"]}"')
		BashSubprocess(Name = f"BWA Index", Command = f'bwa index "{FNdict["FASTA"]}"')
		BashSubprocess(Name = f"GATK Index", Command = f'{GATK_PATH} CreateSequenceDictionary --VERBOSITY ERROR -R "{FNdict["FASTA"]}"')
		BashSubprocess(Name = f'Genome BED', Command = f'awk \'BEGIN {{ FS = "\\t" }}; {{ print $1 FS "0" FS $2 }}\' "{FaiFN}" > "{BedFN}"')
		JsonFN = os.path.join(TempDir, f'{GenomeName}.info.json')
		ConfigJson = CreateGenomeInfo(GenomeName, OutputDir)
		SaveJSON(ConfigJson, JsonFN)
		BashSubprocess(Name = f"Copy Files", Command = f'cp {os.path.join(TempDir, "*")} "{OutputDir}"')

## ------======| CAPTURE PREP |======------

def CreateCaptureInfo(CaptureName, Dir):
	ConfigJson = {
		'NAME':     CaptureName,
		'CAP':      f'{CaptureName}.capture.bed',
		'NOTCAP':   f'{CaptureName}.not_capture.bed'
	}
	ConfigJson = { Key: os.path.join(Dir, Value) for Key, Value in ConfigJson }
	return ConfigJson

def CapturePreparation(CaptureName, InputBED, GenomeInfoJSON, ParentDir):
	OutputDir = os.path.join(ParentDir, CaptureName)
	os.mkdir(OutputDir)
	logging.info(f'Directory created: "{OutputDir}"')
	GenomeInfo = LoadJSON(GenomeInfoJSON)
	with tempfile.TemporaryDirectory() as TempDir:
		FNdict = CreateCaptureInfo(CaptureName, TempDir)
		BashSubprocess(Name = f"CapturePreparation.Filter&Sort", Command = f'set -o pipefail; bedtools intersect -a "{GenomeInfo["BED"]}" -b "{InputBED}" | bedtools sort -faidx "{GenomeInfo["FAI"]}" | sed -e \'s/$/\\t\\./\' > "{FNdict["CAP"]}"')
		BashSubprocess(Name = f"CapturePreparation.NotCapture", Command = f'bedtools subtract -a "{GenomeInfo["BED"]}" -b "{FNdict["CAP"]}" | sed -e \'s/$/\\t\\./\' > "{FNdict["NOTCAP"]}"')
		JsonFN = os.path.join(TempDir, f'{CaptureName}.info.json')
		ConfigJson = CreateCaptureInfo(CaptureName, OutputDir)
		SaveJSON(ConfigJson, JsonFN)
		BashSubprocess(Name = f"Copy Files", Command = f'cp {os.path.join(TempDir, "*")} "{OutputDir}"')
