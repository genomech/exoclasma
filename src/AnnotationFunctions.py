from SharedFunctions import *

def ANNOVAR(InputVCF, OutputVCF, AnnovarFolder, GenomeAssembly, Threads = multiprocessing.cpu_count()):
	logging.info(f'Input VCF: "{InputVCF}"; Output VCF: "{OutputVCF}"; Genome Assembly: {GenomeAssembly}') 
	with tempfile.TemporaryDirectory() as TempDir:
		TableAnnovarPath = os.path.join(AnnovarFolder, "table_annovar.pl")
		AnnovarDBPath = os.path.join(AnnovarFolder, "humandb")
		TempVCF = os.path.join(TempDir, "temp.vcf")
		AnnotatedVCF = f"{TempVCF}.{GenomeAssembly}_multianno.vcf"
		BashSubprocess(Name = f"ANNOVAR.TempVCF", Command = f'cp "{InputVCF}" "{TempVCF}"')
		BashSubprocess(
			Name = f"ANNOVAR.Annotation",
			Command = f'perl "{TableAnnovarPath}" "{TempVCF}" "{AnnovarDBPath}" --buildver {GenomeAssembly} --protocol refGene --operation g --remove --vcfinput --thread {Threads}',
			AllowedExitCodes = [25]
		)
		BashSubprocess(Name = f"ANNOVAR.BgzipVCF", Command = f'bgzip -c "{AnnotatedVCF}" > "{OutputVCF}"')
		BashSubprocess(Name = f"ANNOVAR.IndexVCF", Command = f'tabix -p vcf "{OutputVCF}"')
