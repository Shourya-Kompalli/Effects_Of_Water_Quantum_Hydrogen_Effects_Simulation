#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
command -v i-pi >/dev/null || { echo "i-pi is not installed in this environment." >&2; exit 1; }
command -v lmp >/dev/null || { echo "LAMMPS executable 'lmp' is not installed or not on PATH." >&2; exit 1; }

python data/prepare_lammps_data.py

i-pi data/input.xml > simulation.log 2>&1 &
ipi_pid=$!
worker_pids=()

cleanup() {
    for pid in "${worker_pids[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
            kill "$pid" 2>/dev/null || true
        fi
    done
    if kill -0 "$ipi_pid" 2>/dev/null; then
        kill "$ipi_pid" 2>/dev/null || true
    fi
}
trap cleanup EXIT INT TERM

sleep 2
if ! kill -0 "$ipi_pid" 2>/dev/null; then
    cat simulation.log >&2
    echo "i-PI exited before the LAMMPS clients connected." >&2
    exit 1
fi

for worker in 1 2; do
    lmp -log "lammps-${worker}.log" -in data/in.lmp &
    worker_pids+=("$!")
done

wait "$ipi_pid"
for pid in "${worker_pids[@]}"; do
    wait "$pid"
done
