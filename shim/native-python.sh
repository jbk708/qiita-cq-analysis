#!/bin/bash
# SLURM_NATIVE_PYTHON wrapper: points ingest_ena_reads jobs at the ENA shim while portal.env exists.
# The path must stay space-free (the orchestrator interpolates it unquoted).
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
STATE=${SHIM_STATE_DIR:-/ddn_scratch/$(id -un)/qiita-pilot/ena-shim}
PY=${SHIM_PY:-/ddn_scratch/$(id -un)/qiita-dev/Qiita/qiita-compute-orchestrator/.venv/bin/python}

if [[ " $* " == *" ingest_ena_reads "* && -f $STATE/portal.env ]]; then
    verdict=$("$PY" -I "$HERE/scope_check.py" "$STATE/portal.runs") || {
        echo "ena-shim: scope check failed" >&2
        exit 1
    }
    # A proxy-mode shim forwards the runs it doesn't hold to ENA, so a mixed job can use it.
    if [[ $verdict == mixed ]] && grep -qx ENA_SHIM_MODE=proxy "$STATE/portal.env"; then
        verdict=in
    fi
    case $verdict in
        in)
            set -a
            . "$STATE/portal.env"
            set +a
            [[ ${MIINT_ENA_PORTAL_URL_BASE:-} == http://* ]] || {
                echo "ena-shim: portal.env has no http:// MIINT_ENA_PORTAL_URL_BASE" >&2
                exit 1
            }
            export QIITA_DEV_ENA_PRIORITY_RUNS=$STATE/portal.runs
            echo "ena-shim: portal=$MIINT_ENA_PORTAL_URL_BASE" >&2
            ;;
        out) ;;
        *)
            echo "ena-shim: job mixes shim and non-shim runs ($verdict)" >&2
            exit 1
            ;;
    esac
fi
exec "$PY" "$@"
