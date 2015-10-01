#-----------------------------------------------------------------------------
#
# By (un)commenting lines below, can choose between Intel + MKL or GNU + gsl
# Be sure to not mix up those up; changing the compiler requires changing the
# libs and most of the optimizations.
#
# Started: Fri 07 Sep 2012, 01:32 PM
# Modified: Tue 09 Oct 2012, 02:34 PM
# Modified: Mon 07 Jan 2013, 02:12 PM spelling correction
# Last Modified: Mon 04 Nov 2013, 06:13 AM
#-----------------------------------------------------------------------------

# F90 = ifort
F90 = gfortran
CC = gcc

F95 = $(F90)
OPTS = -O3 

# MKLROOT = /opt/intel/composer_xe_2013.3.163/mkl
LIBS =   -L$(MKLROOT)/lib/intel64 -lmkl_rt -lpthread -lm
INCS = -I.

#---------------------------------------
# For GNU gfortran and GSL: 
#---------------------------------------
# LIBS = -L/usr/lib/ -lgsl -lblas
#-----------------------------------------------------------------------------
#
CC = gcc

#LIBS = -L/usr/lib/libblas/ -lblas

