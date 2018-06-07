# Build settings for Cray compiler with OpenACC
export F90=ftn

# f2py does not break long lines so set the maximum
# line length to the maximum permitted by the Cray compiler
F90FLAGS="-N 1023"

#F90FLAGS+="-g -O0 -h noomp -Rb -K trap=fp,denorm"
F90FLAGS+="-O2"
#F90FLAGS+="-O3 -O ipa5 -rd"
#-h loop_trips=medium

# Turn on OpenACC
# Set $CRAY_ACC_DEBUG to 1, 2 or 3 to get run-time debug output
# of accelerator activity
F90FLAGS+="-h acc"

# Display compiler options in effect and generate annotated
# listing
F90FLAGS+="-h display_opt -rm"
#F90FLAGS+="-h profile_generate"

LDFLAGS=
#LDFLAGS="-h profile_generate"

# Flags for whole-program optimisation
#F90FLAGS+="-hpl=${PWD}/nemolite2d_wpl_dir -h wp"
#LDFLAGS+="-hpl=${PWD}/nemolite2d_wpl_dir"

AR=ar

export F90FLAGS
export LDFLAGS
export AR
