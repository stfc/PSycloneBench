# Build settings for Gnu compiler
F90 = gfortran
CC = gcc

NETCDF_INC = /usr/include
NETCDF_LIB = /usr/lib64

#CFLAGS = -g -O0
CFLAGS = -O3 -Wall -Wextra -std=c11

CFLAGS +=  -I${NETCDF_INC}

LDFLAGS += -L${NETCDF_LIB} -lm

