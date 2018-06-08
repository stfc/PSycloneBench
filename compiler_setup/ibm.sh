# Build settings for ibm compiler
F90=xlf2003_r
CC=xlc

CFLAGS="-O3"

# Ideally would specify -qlanglvl=2003std but GungHo/src/utils/utils.F90
# won't build with that option :-(

F90FLAGS="-qcheck -g -O0"
#F90FLAGS="-O3"
LDFLAGS=

AR=ar

export F90
export F90FLAGS
export CC
export CFLAGS
export LDFLAGS
export AR
