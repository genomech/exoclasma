from Bio import SeqIO
import argparse
import bz2
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
import warnings

## ------======| LOGGING |======------

def ExceptionHook(Type, Value, Traceback): logging.exception(f"{Type.__name__}: {Value}", exc_info=(Type, Value, Traceback))

def ConfigureLogger(LogFileName=os.devnull, Level=logging.DEBUG):
	Formatter = "%(asctime)-30s%(levelname)-13s%(funcName)-35s%(message)s"
	logger = logging.getLogger()
	while logger.hasHandlers(): logger.removeHandler(logger.handlers[0])
	Handlers = [logging.FileHandler(LogFileName), logging.StreamHandler(sys.stderr)]
	logging.basicConfig(level=Level, format=Formatter, handlers=Handlers)
	sys.excepthook = ExceptionHook
	
def Timestamp(TS): return datetime.timedelta(seconds=(time.time() - TS))


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


## ------======| THREADING |======------

@contextlib.contextmanager
def Threading(Name, Threads = multiprocessing.cpu_count()):
	StartTime = time.time()
	pool = multiprocessing.Pool(Threads)
	yield pool
	pool.close()
	pool.join()
	del pool
	logging.info(f"Name: {Name}; threads: {Threads}; total time: {Timestamp(StartTime)}")


## ------======| SUBPROCESS |======------

def BashSubprocess(Name, Command, AllowedExitCodes = list()):
	StartTime = time.time()
	logging.debug(Command)
	Shell = subprocess.Popen(Command, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	Stdout, Stderr = Shell.communicate()
	if (Shell.returncode != 0) and (Shell.returncode not in AllowedExitCodes):
		logging.error(f"Name: {Name}; exit code: {Shell.returncode}; shell command: '{Command}'; details:\n{Stderr.decode('utf-8')}")
		raise RuntimeError
	if Shell.returncode in AllowedExitCodes: logging.warning(f"Name: {Name}; exit code: {Shell.returncode} [ALLOWED]")
	logging.info(f"Name: {Name}; total time: {Timestamp(StartTime)}")
	return Stdout


## ------======| LOAD CONFIG |======------

CONFIG_RESTRICTION_ENZYMES = LoadJSON(os.path.join(CurrentDir(), "..", "config", "RestrictionEnzymes.json"))
CONFIG_JAVA_OPTIONS = LoadJSON(os.path.join(CurrentDir(), "..", "config", "JavaOptions.json"))
GATK_PATH = os.path.abspath(os.path.join(CurrentDir(), "..", "gatk", "gatk"))
