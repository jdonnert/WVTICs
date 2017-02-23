SHELL = /bin/bash

## OPTIONS  ##
#OPT     += -DSAVE_WVT_STEPS         # write IC file for every WVT step
#OPT     += -DSPH_CUBIC_SPLINE      # for use with Gadget2
OPT     += -DREJECTION_SAMPLING      # use von Neumann rejection sampling to improve initial random positions


ifndef SYSTYPE
    SYSTYPE := $(shell hostname)
endif

OPTIMIZE = -Wall -O3
GSL_INCL = $(CPPFLAGS)
GSL_LIBS = $(LD_LIBRARY_FLAGS)

## TARGET ##

EXEC = WVTICs

## FILES ##

SRCDIR    = src/

SRCFILES := ${shell find $(SRCDIR) -name \*.c -print} # all .c files in SRCDIR
OBJFILES = $(SRCFILES:.c=.o)

INCLFILES := ${shell find src -name \*.h -print} # all .h files in SRCDIR
INCLFILES += Makefile

CFLAGS  = -std=c99 -fopenmp $(OPTIMIZE) $(OPT) $(GSL_INCL) $(FFTW_INCL)

LINK    = $(GSL_LIBS) -lm -lgsl -lgslcblas 

## RULES ##

%.o : %.c
	@echo [CC] $@
	@$(CC) $(CFLAGS)  -o $@ -c $<

$(EXEC) : $(OBJFILES)
	@echo $(CC) 
	$(CC) $(CFLAGS) $(OBJFILES) $(LINK) -o $(EXEC)
	@ctags -w $(SRCFILES) $(INCLFILES)

$(OBJFILES) : $(INCLFILES) $(SRCFILES)

clean : 
	rm $(OBJFILES) $(EXEC)
