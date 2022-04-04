#!/bin/bash

Exoclasma=$(realpath ../exoclasma)
BED=$(realpath ./testdata/TestCapture_sacCer3_Exome.bed)
Name=".test.CapturePreparation"
GenomeInfo=$(realpath ./.test.RefseqPreparation/.test.RefseqPreparation.info.json)
Parent=$(realpath ./)
TargetDir=""$Parent"/"$Name""

rm -r "$TargetDir"

python3 "$Exoclasma" CapturePreparation --bed "$BED" --name "$Name" --parent "$Parent" --genomeinfo "$GenomeInfo"
