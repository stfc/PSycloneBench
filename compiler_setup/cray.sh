# Build settings for cray compiler
F90=ftn
CC=cc

CFLAGS="-O3 -h c99 -h aggress"

# f2py does not break long lines so set the maximum
# line length to the maximum permitted by the Cray compiler
F90FLAGS="-N 1023"

#F90FLAGS+=" -g -O0 -h noomp -Rb -K trap=fp,denorm"
F90FLAGS+=" -O3 -O ipa5 -rd"
#-h loop_trips=medium
# Have to explicitly turn OpenMP off if building 
# with optimising flags
F90FLAGS+=" -h noomp"

# Display compiler options in effect and generate annotated
# listing
F90FLAGS+=" -h display_opt -rm"
F90FLAGS+=" -h profile_generate"

LDFLAGS="-h profile_generate"

# Flags for whole-program optimisation
F90FLAGS+=" -hpl=${PWD}/nemolite2d_wpl_dir -h wp"
LDFLAGS+=" -hpl=${PWD}/nemolite2d_wpl_dir"

AR=ar

export F90
export F90FLAGS
export CC
export CFLAGS
export LDFLAGS
export AR
