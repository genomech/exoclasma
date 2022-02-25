from src.SharedFunctions import *

def SeriesReplace(String, Replacements):
	Rep = dict((re.escape(key), value) for key, value in Replacements.items()) 
	return re.compile("|".join(Rep.keys())).sub(lambda x: Rep[re.escape(x.group(0))], String)

def RefseqPreparation(FastaPath, GenomeName, ParentDir):
	OutputDir = os.path.join(ParentDir, GenomeName)
	os.mkdir(OutputDir)
	logging.info(f"Directory created: \"{OutputDir}\"")
	Fasta = SeqIO.parse(OpenAnyway(FastaPath, "rt"), "fasta")
	logging.info(f"FASTA opened: \"{FastaPath}\"")
	SearchQueries = { name: re.compile("|".join([f"({SeriesReplace(item, CONFIG_RESTRICTION_ENZYMES['WildCards'])})" for item in sequences])) for name, sequences in CONFIG_RESTRICTION_ENZYMES["Enzymes"].items() }
	with tempfile.TemporaryDirectory() as TempDir:
		NewFastaFN = os.path.join(TempDir, f"{GenomeName}.fasta")
		ChromSizesFN = os.path.join(TempDir, f"{GenomeName}.chromsizes.txt")
		with open(NewFastaFN, "w") as NewFasta, open(ChromSizesFN, "w") as ChromSizes:
			RestrictionSites = { name: os.path.join(TempDir, f"{GenomeName}.restriction_sites.{name}.txt") for name, sequences in CONFIG_RESTRICTION_ENZYMES["Enzymes"].items() }
			for contig in Fasta:
				contig.name = re.sub("[^\w\.]", "_", contig.name)
				Seq = contig.seq.__str__()
				SeqLength = len(Seq)
				SeqIO.write([contig], NewFasta, "fasta")
				ChromSizes.write(f"{contig.name}\t{SeqLength}\n")
				for enzyme, query in SearchQueries.items():
					FileWrapper = open(RestrictionSites[enzyme], "a")
					FileWrapper.write(" ".join([contig.name] + [str(match.start() + 1) for match in query.finditer(Seq)] + [str(SeqLength)]) + "\n")
					FileWrapper.close()
				logging.debug(f"Contig \"{contig.name}\" ready")
		logging.info(f"New FASTA, ChromSizes and RestrictionSites are ready")
		BashSubprocess(Name = f"SAMtools Index", Command = f"samtools faidx \"{NewFasta.name}\"")
		BashSubprocess(Name = f"BWA Index", Command = f"bwa index \"{NewFasta.name}\"")
		BashSubprocess(Name = f"GATK Index", Command = f"gatk/gatk CreateSequenceDictionary -R \"{NewFasta.name}\"")
		BashSubprocess(Name = f"Copy Files", Command = f"cp {os.path.join(TempDir, '*')} \"{OutputDir}\"")
