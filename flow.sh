#!/bin/sh
conda activate fipy34;

mpirun -np 8 python finitevolume.py;

foldername=$(date +%Y%m%d_%H%M%S);
mkdir "$foldername";
cp params.yml ./"$foldername";
mv *.vtk ./"$foldername";
mv nohup* ./"$foldername"
