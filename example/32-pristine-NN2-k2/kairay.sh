#!/bin/bash
#SBATCH -J RW
#SBATCH -o task.%j.out
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -p debug

module purge
module add intel/17.0.5.239 impi/2017

CASE="bivo4"
for i in 3200 12800
do
    DIR=${CASE}-t$i
    if [ ! -d $DIR ] 
    then
        sed -i "3s/.*/$i/" ${CASE}.in
        /home/wufeng/work/BiVO4Mo-RW/modelvssm/build/calc_activation.py $CASE -o $DIR -n 16 --np 16
    fi
done
