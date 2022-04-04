#!/bin/bash

TargetDir=".test.AlignPipeline"
ExoCLasma="../exoclasma"
Units="Test.AlignPipeline.Units.json"

rm -r ""$(realpath $TargetDir)""
python3 ""$(realpath $ExoCLasma)"" AlignPipeline -u ""$(realpath $Units)""
