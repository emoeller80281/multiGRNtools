#!/bin/bash
# Single source of truth for how many resources an inference method may use.
#
# The benchmark only means something if every method is handed the same budget, so that
# budget is defined in exactly one place -- the #SBATCH header of the top-level runner
# (run_tools.sh / run_tools_stability.sh) -- and read back out of the environment here.
# Method scripts must not decide for themselves: they consume ALLOC_CPUS and ALLOC_MEM_*.
#
# Before this existed each method picked its own number, so the "comparison" was between
# methods running on different hardware budgets. The worst case was CellOracle, which
# hard-coded n_jobs=4 while every other method took the full allocation.
#
# Exports:
#   ALLOC_CPUS      cores this run may use          (from --cpus-per-task / -c)
#   ALLOC_MEM_MB    memory this run may use, MB     (from --mem / --mem-per-cpu)
#   ALLOC_MEM_GB    the same figure in whole GB, for tools that want GB
#
# Plus the thread-pool variables below, so libraries that never see ALLOC_CPUS still stay
# inside the budget.
#
# Sourced, not executed:  source "${PROJECT_DIR}/src/common/resource_env.sh"

# ===== CORES =====
# SLURM_CPUS_PER_TASK is what -c / --cpus-per-task sets. Outside Slurm, `nproc` reports the
# cgroup/affinity-limited count, which is what we want -- note `nproc --all` does NOT: on
# these nodes it reports all 88 physical cores regardless of the allocation, which is exactly
# the over-subscription this file exists to prevent.
if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
    ALLOC_CPUS="${SLURM_CPUS_PER_TASK}"
elif command -v nproc > /dev/null 2>&1; then
    ALLOC_CPUS="$(nproc)"
else
    ALLOC_CPUS=1
fi

# ===== MEMORY =====
# --mem sets SLURM_MEM_PER_NODE (MB); --mem-per-cpu sets SLURM_MEM_PER_CPU (MB) instead, so
# multiply that one back up by the core count to get the real per-task budget.
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    ALLOC_MEM_MB="${SLURM_MEM_PER_NODE}"
elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
    ALLOC_MEM_MB=$(( SLURM_MEM_PER_CPU * ALLOC_CPUS ))
elif [ -r /proc/meminfo ]; then
    # Not under Slurm: fall back to physical RAM so the variables are always populated.
    ALLOC_MEM_MB=$(awk '/^MemTotal:/ { printf "%d", $2 / 1024 }' /proc/meminfo)
else
    ALLOC_MEM_MB=0
fi
ALLOC_MEM_GB=$(( ALLOC_MEM_MB / 1024 ))

export ALLOC_CPUS ALLOC_MEM_MB ALLOC_MEM_GB

# ===== IMPLICIT THREAD POOLS =====
# Several methods never pass a thread count to their numeric libraries, which then size their
# own pools. OpenBLAS happens to honour the cpuset here (measured: 12 threads inside a 12-core
# allocation), but that is a property of this build rather than a guarantee, and MKL and
# OpenMP elsewhere size from the physical core count. Pinning them makes the budget explicit
# and identical for every method instead of depending on which BLAS an env happens to ship.
#
# ALLOC_BLAS_THREADS exists because there is a real trade-off here and no single right answer.
# A method that forks N workers which each open an ALLOC_CPUS-wide BLAS pool oversubscribes
# (the cpuset still caps real parallelism at ALLOC_CPUS, so this costs that method time rather
# than stealing from anyone else). Set ALLOC_BLAS_THREADS=1 before sourcing to use the usual
# HPC convention of one BLAS thread per forked worker.
ALLOC_BLAS_THREADS="${ALLOC_BLAS_THREADS:-${ALLOC_CPUS}}"

export ALLOC_BLAS_THREADS
export OMP_NUM_THREADS="${ALLOC_BLAS_THREADS}"
export OPENBLAS_NUM_THREADS="${ALLOC_BLAS_THREADS}"
export GOTO_NUM_THREADS="${ALLOC_BLAS_THREADS}"
export MKL_NUM_THREADS="${ALLOC_BLAS_THREADS}"
export NUMEXPR_NUM_THREADS="${ALLOC_BLAS_THREADS}"
export VECLIB_MAXIMUM_THREADS="${ALLOC_BLAS_THREADS}"
# data.table and R's parallel::mclapply read these two respectively; R's mc.cores option is
# initialised from MC_CORES, which is what DIRECTNET's mclapply backend ends up using.
export R_DATATABLE_NUM_THREADS="${ALLOC_BLAS_THREADS}"
export MC_CORES="${ALLOC_CPUS}"

# ===== REPORT =====
# Printed once per run so the .out file records the budget the numbers were produced under.
resource_env_summary() {
    local source_desc
    if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
        source_desc="#SBATCH -c ${SLURM_CPUS_PER_TASK} in the submitted runner"
    else
        source_desc="not under Slurm, detected from nproc"
    fi

    echo "===== RESOURCE ALLOCATION ====="
    echo "  cores        : ${ALLOC_CPUS}"
    echo "  memory       : ${ALLOC_MEM_MB} MB (${ALLOC_MEM_GB} GB)"
    echo "  blas threads : ${ALLOC_BLAS_THREADS}"
    echo "  source       : ${source_desc}"
}
