from Bio import SeqIO
import argparse
import bisect
import bz2
import collections
import contextlib
import copy
import datetime
import functools
import glob
import gzip
import io
import json
import logging
import math
import multiprocessing
import numpy
import os
import pandarallel
import pandas
import pysam
import re
import subprocess
import sys
import tabix
import tempfile
import time
import types
import vcf
import warnings


## ------======| LOGGING |======------

VERBOSITY_LEVELS = {'debug': 10, 'info': 20, 'warning': 30, 'error': 40}

def ExceptionHook(Type, Value, Traceback): logging.exception(f"{Type.__name__}: {Value}", exc_info=(Type, Value, Traceback))

def ConfigureLogger(LogFileName=os.devnull, Level=logging.DEBUG):
	Formatter = "%(asctime)-30s%(levelname)-13s%(funcName)-35s%(message)s"
	logger = logging.getLogger()
	while logger.hasHandlers(): logger.removeHandler(logger.handlers[0])
	Handlers = [logging.FileHandler(LogFileName), logging.StreamHandler(sys.stderr)]
	logging.basicConfig(level=Level, format=Formatter, handlers=Handlers)
	sys.excepthook = ExceptionHook
	
def Timestamp(TS): return str(datetime.timedelta(seconds=(time.time() - TS)))[:-7]

def RenderParameters(Parameters):
	for Key, Value in Parameters.items():
		if type(Value) == list: Value = '; '.join(Value)
		logging.info(f'{Key.replace("_", " ")}: {Value}')
	return types.SimpleNamespace(**Parameters)

## ------======| I/O |======------

def CurrentDir(): return os.path.dirname(os.path.abspath(__file__))

def GzipCheck(FileName): return open(FileName, 'rb').read(2).hex() == "1f8b"

def Bzip2Check(FileName): return open(FileName, 'rb').read(3).hex() == "425a68"

def CompressedOpenFunc(Type): return {"none": open, "gzip": gzip.open, "bz2": bz2.open}[Type]

def OpenAnyway(FileName, Mode):
	CheckFlags = GzipCheck(FileName = FileName), Bzip2Check(FileName = FileName)
	Type = {(1, 0): "gzip", (0, 1): "bz2", (0, 0): "none"}[CheckFlags]
	return CompressedOpenFunc(Type)(FileName, Mode)


## ------======| JSON |======------

def SaveJSON(Data, FileName, Compression = "none"):
	return json.dump(Data, CompressedOpenFunc(Compression)(FileName, 'wt'), indent = 4, ensure_ascii = False)

def LoadJSON(FileName): return json.load(OpenAnyway(FileName, 'rt'))

def RecursiveDict(PathDict):
	Result = {}
	Flat = [key for key in PathDict.keys() if len(key) == 1]
	for key in Flat: Result[key[0]] = PathDict[key]
	Unique = list(set([key[0] for key in PathDict.keys() if len(key) > 1]))
	for Key in Unique: Result[Key] = RecursiveDict({tuple(key[1:]): value for key, value in PathDict.items() if key[0] == Key})
	return Result

## ------======| THREADING |======------

@contextlib.contextmanager
def Threading(Name, Threads = multiprocessing.cpu_count()):
	logging.info('* entry point *')
	StartTime = time.time()
	pool = multiprocessing.Pool(Threads)
	yield pool
	pool.close()
	pool.join()
	del pool
	logging.info(f'{Name} finished on {Threads} threads - {Timestamp(StartTime)}')


## ------======| SUBPROCESS |======------

def BashSubprocess(Name, Command, AllowedExitCodes = list()):
	StartTime = time.time()
	logging.debug(f'CMD: {Command}')
	Shell = subprocess.Popen(Command, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	Stdout, Stderr = Shell.communicate()
	if (Shell.returncode != 0) and (Shell.returncode not in AllowedExitCodes):
		for Line in [f'{Name} finished with non-zero exit code: {Shell.returncode}', f'Shell command: "{Command}"', f'Details:\n{Stderr.decode("utf-8")}']: logging.error(Line)
		raise RuntimeError
	if Shell.returncode in AllowedExitCodes: logging.warning(f'{Name} finished with ALLOWED non-zero exit code: {Shell.returncode}')
	logging.info(f'{Name} finished - {Timestamp(StartTime)}')
	return Stdout

## ------======| LOAD CONFIG |======------

CONFIG_RESTRICTION_ENZYMES = LoadJSON(os.path.join(CurrentDir(), "..", "config", "RestrictionEnzymes.json"))

CONFIG_JAVA_OPTIONS = LoadJSON(os.path.join(CurrentDir(), "..", "config", "JavaOptions.json"))

CONFIG_ADAPTERS = LoadJSON(os.path.join(CurrentDir(), "..", "config", "Adapters.json"))

GATK_PATH = os.path.abspath(os.path.join(CurrentDir(), "..", "gatk", "gatk"))

JUICERTOOLS_PATH = "/Data/Tools/juicer_tools_1.22.01.jar" # TODO

GLOBAL_BACKUP = None
