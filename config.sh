# Paths for the pilot scripts on barnacle. Override any of these in the environment.
export PILOT_REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export QDEV_ROOT=${QDEV_ROOT:-/ddn_scratch/$USER/qiita-dev}
export QDEV_SHARE=${QDEV_SHARE:-/ddn_scratch/$USER/qiita-dev-share}
export QDEV_CONDA=${QDEV_CONDA:-/ddn_scratch/$USER/conda-envs/qiita-reads}
export QDEV_ACCOUNT=${QDEV_ACCOUNT:-}
# qiita-dev.sh from jbk708/qiita-barnacle-dev
export QDEV_CLI=${QDEV_CLI:-/ddn_scratch/$USER/qiita-barnacle-dev/stack/qiita-dev.sh}
