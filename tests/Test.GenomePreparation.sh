#!/bin/bash

Exoclasma=$(realpath ../exoclasma)
Fasta=$(realpath ./testdata/TestGenome_sacCer3.fa.gz)
Name=".test.RefseqPreparation"
Parent=$(realpath ./)
TargetDir=""$Parent"/"$Name""

rm -r "$TargetDir"

python3 "$Exoclasma" RefseqPreparation --fasta "$Fasta" --name "$Name" --parent "$Parent"
