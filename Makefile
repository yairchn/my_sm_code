# Makefile for various platforms
# Execute using Build csh-script only!
# Used together with Perl scripts in SRC/SCRIPT 
# (C) 2005 Marat Khairoutdinov
#------------------------------------------------------------------

##INC_NETCDF   := /home/software/eaps/pgi/15.10/pkg/netcdf/3.6.2/include/
##LIB_NETCDF   := /home/software/eaps/pgi/15.10/pkg/netcdf/3.6.2/lib
SAM = SAM_$(RAD_DIR)

# Determine platform 
PLATFORM := $(shell uname -s)


# Linux, Yellowstone, Cheyenne
# #

 ifeq ($(PLATFORM),Linux)

#FF77 = mpipf90 -c -Mextend
#FF90 = mpipf90 -c -Mfreeform
#CC = mpipcc -c  -DLINUX
#
#FFLAGS = -Mnoframe -Mvect -Munroll -O2 -Mbyteswapio  
#
#FFLAGS += -I${INC_NETCDF}
#LD = mpipf90
#LDFLAGS =  -L${LIB_NETCDF} -lnetcdf


# FF77 = ftn -fpic -i_dynamic -mcmodel=large -c -fixed -extend_source
# FF90 = ftn -fpic -i_dynamic -mcmodel=large -c -free
# CC = mpcc -fpic -i_dynamic -mcmodel=large -c -O2 -DLINUX
 FF77 = mpif90 -c -fixed -extend_source
 FF90 = mpif90 -c -free
 CC = mpicc -c -O3 -DLINUX


# FFLAGS = -xHost -O3 -pad 
FFLAGS =  -O2 
#FFLAGS =  -O1  -mcmodel=large  
#FFLAGS =  -O2  -xCORE-AVX2  -mcmodel=large  
# FFLAGS = -g -O3 -fpe0 -traceback -pad  -mcmodel=large
# FFLAGS = -g -ftrapuv -check all -traceback -W1 

## FFLAGS += -I${INC_NETCDF}
# LD = mpif90 -fpic -i_dynamic -mcmodel=large
 LD = mpif90  -mcmodel=large
## LDFLAGS = -L${LIB_NETCDF} -lnetcdff 


 endif



# you dont need to edit below this line


#compute the search path
dirs := . $(shell cat Filepath)
VPATH    := $(foreach dir,$(dirs),$(wildcard $(dir))) 

.SUFFIXES:
.SUFFIXES: .f .f90 .c .o



all: $(SAM_DIR)/$(SAM)


SOURCES   := $(shell cat Srcfiles)

Depends: Srcfiles Filepath
	$(SAM_SRC)/SCRIPT/mkDepends Filepath Srcfiles > $@

Srcfiles: Filepath
	$(SAM_SRC)/SCRIPT/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES))) 

$(SAM_DIR)/$(SAM): $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)


.f90.o:
	${FF90}  ${FFLAGS} $<
.f.o:
	${FF77}  ${FFLAGS} $<
.c.o:
	${CC}  ${CFLAGS} -I$(SAM_SRC)/TIMING $(NOTIMERS) $<



include Depends



clean: 
	rm ./OBJ/*



