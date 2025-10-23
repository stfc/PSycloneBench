# Build settings for the Nvidia compiler
# ================================================
# Fortran compiler

# ==============================
export F90=$FC
export PSYCLONE_NVIDIA_LIB_DIR=${HOME}/Projects/PSyclone/lib/profiling/nvidia
export OMPTARGETFLAGS="-mp=gpu -gpu=ccnative"
export OMPFLAGS="-mp"
export UMEMFLAGS="-gpu=managed"

