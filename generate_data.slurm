#!/bin/bash
#SBATCH --job-name="accML seed generation"
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --partition=gpu

echo "Starting"
module load miniforge3
module load matlab
conda activate torchrl
echo "Modules loaded, starting script"
cd /home/lqiao/vf38/acceleratorML/
python generate_data.py
echo "Job finished"