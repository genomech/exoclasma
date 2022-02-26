from src.SharedFunctions import *

def VcfBlacklist(VcfList, OutputVCF, Threads = multiprocessing.cpu_count()):
	VCFLog = ', '.join([f'"{item}"' for item in VcfList])
	logging.info(f'VCF List: "{VCFLog}"; Output VCF: "{OutputVCF}"')
