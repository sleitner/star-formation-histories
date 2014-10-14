MYGSL_DIR=/home/sleitner/simplemsi/include/

F90FLAGS = -O
CC = icc
CFLAGS = -O3 -mcmodel=medium -i_dynamic -w1
CFLAGS = -mcmodel=medium -i_dynamic -w1
F90 = gfortran
CF90 = gfortran

#need gsl
LIBS = -lm -lgsl -lgslcblas # -L /opt/intel/Compiler/11.0/081/lib/intel64/ -lirc -lcxaguard
INCLUDES = -I$(MYGSL_DIR) -L$(MYGSL_DIR)/lib/
PROGS = simple_msi.x

DOTH = hfiles/auxiliary.h hfiles/cosmology.h hfiles/defs.h
OBJS = auxiliary.c cosmology.c 
CHECKFILES = $(OBJS) $(DOTH)

all : $(PROGS)

#-----------------------------------------------
simple_msi.x: $(CHECKFILES) simple_msi.c
	$(CF90) -I. $(INCLUDES) $(OBJS) simple_msi.c \
		-o simple_msi.x \
		$(LIBS)

#-----------------------------------------------

clean: 
	rm -f $(PROGS)


