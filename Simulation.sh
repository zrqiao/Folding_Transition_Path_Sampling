#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-10:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shaknovich       # Partition to submit to
#SBATCH --mem=1000          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o kmcoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e kmcerrors_%j.err  # File to which STDERR will be written, %j inserts jobid
export FEL_NAME=$1

cp simulate_tp_1d.py $FEL_NAME\/
cd $FEL_NAME\/
python ./simulate_tp_1d.py PES.dat --stored-paths tps_$FEL_NAME\.dat
cd ../
python ./Transit Time Plot.py
python ./Displacement_Plot.py

