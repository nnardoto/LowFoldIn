#!/bin/bash
#SBATCH -p ep-short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 01:30:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
# A25: termos de face (a25_face.py) no espelho gen/: calib -> x1 -> check -> tr
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export A25_ROOT=$SLURM_SUBMIT_DIR/gen/
PY=/home/wendel/packages/miniconda3/bin/python3; S=$A25_ROOT/lfi_tese/trincheira_vol/scripts/a25_face.py
cd $SLURM_SUBMIT_DIR
$PY $S calib > calib.txt 2>&1 || exit 1
$PY $S x1 $SLURM_SUBMIT_DIR/face_x1/a25_x1.sitepot > gen_x1.txt 2>&1 || exit 2

$PY $S tr $SLURM_SUBMIT_DIR/trench/a25_tr.sitepot > gen_tr.txt 2>&1 || exit 3
echo fim
