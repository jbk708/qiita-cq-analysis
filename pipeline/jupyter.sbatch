#!/bin/bash
# JupyterLab on barnacle's jupyter partition with the qiita-reads env.
# Submit: sbatch pipeline/jupyter.sbatch (QDEV_CONDA overrides the env)
# Then read the log for the node, port and token, and tunnel from your laptop:
#   ssh -N -L 8888:<node>:<port> barnacle      ->  open http://localhost:8888/?token=<token>
#SBATCH -J jupyter-qiita-reads
#SBATCH -p jupyter
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -o jupyter-%j.log
ENV=${QDEV_CONDA:-/ddn_scratch/$USER/qiita-pilot/conda-envs/qiita-reads}
PORT=$(shuf -i 20000-29999 -n 1)
echo "node: $(hostname)  port: $PORT"
echo "tunnel: ssh -N -L 8888:$(hostname):$PORT barnacle"
exec $ENV/bin/jupyter lab --no-browser --ip=0.0.0.0 --port=$PORT --notebook-dir=$HOME
