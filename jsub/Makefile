MAKEFILE      = Makefile
# local variables
OPTFLAGS      = -O4 -O5 -ffast-math
WARNFLAGS     = -Wshadow
DBGFLAGS      = -g3
XTRAFLAGS     = -fsigned-char -Wno-free-nonheap-object
# global variables
AR	      = ar
COMPILER      = g++
INCLUDE       = -I../cmlib -I../include -I../atmos -I../grid
COMPFLAGS     = $(OPTFLAGS) $(WARNFLAGS) $(DBGFLAGS) $(XTRAFLAGS) $(TYPE_DEF)
# export global variables
export MAKE
export COMPILER
export INCLUDE
export COMPFLAGS
# ...
LIBPATH       = -L../cmlib -L../lib -L../grid -L../atmos
LIBS	      = -latmos -lcg -lcm -lz -lm -lpthread
DEST	      = .
PROGRAM       = jsub
SHELL         = /bin/sh

CC	      = $(COMPILER)
CCFLAGS	      = $(COMPFLAGS) $(INCLUDE) $(INCPATH)
LD	      = $(COMPILER)
LDFLAGS	      = $(LIBPATH) $(DBGFLAGS)

OBJS	      = cfg.o jsub.o

.cc.o:
	${CC} -c ${CCFLAGS} $< 

.SUFFIXES: .o .cc

all:
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) all -C $$subdir ; done
	$(MAKE) $(OBJS) $(HEADERS)
	@echo "Linking $(PROGRAM) ..."
	$(LD) -o $(PROGRAM) $(LDFLAGS) $(OBJS) $(LIBS)
	@echo "done"

clean:
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) clean -C $$subdir ; done
	rm -f $(OBJS) core

distclean: 
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) distclean -C $$subdir ; done
	rm -f $(OBJS) core $(PROGRAM) $(MAKEFILE).bak $(MAKEFILE)

install:
	cp -v $(PROGRAM) /usr/local/bin/

dep:
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) dep -C $$subdir ; done
	$(MAKE) makedepend

makedepend:
	@makedepend -f$(MAKEFILE) $(OBJS:.o=.cc) $(INCLUDE)

update:
	$(DEST)/$(PROGRAM) 

# DO NOT DELETE THIS LINE -- make depend depends on it.
