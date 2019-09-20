#!/bin/sh

singularity exec ./ocellaris_2019.0.2.sif mpirun -np 8 python3 binary.py 2>outputerr.txt;

foldername=$(date +%Y%m%d_%H%M%S);
mkdir "$foldername";
cp params.yml ./"$foldername";
mv output* ./"$foldername"
