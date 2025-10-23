# Build settings for the Nvidia compiler
# ================================================
# Fortran compiler

# ==============================
export F90=$FC

export OMPTARGETFLAGS="-mp=gpu -gpu=ccnative"
export OMPFLAGS="-mp"
export UMEMFLAGS="-gpu=mem:managed"
export ACCFLAGS="-acc=gpu -gpu=ccnative"

