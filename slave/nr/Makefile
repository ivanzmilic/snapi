MAKEFILE      = Makefile

#COMPFLAGS     = $(OPTFLAGS) $(WARNFLAGS) $(DBGFLAGS) $(EXTCOMPFLAGS)

INCLUDE	      += -I.. -I../../cmlib -I../../include -I/usr/local/include
DEST	      = .

AR	      = ar
COMPILER      = g++
RANLIB	      = ranlib

CC	      = $(COMPILER)
CCFLAGS	      = $(COMPFLAGS) $(INCLUDE) 

OBJS          = fft_fftw.o
LIBRARY       = libnr.a

.cc.o:
		$(CC) -c $(CCFLAGS) $<

.SUFFIXES:	.o .cc
	
all:            $(LIBRARY)

$(LIBRARY):     $(OBJS) ${LIBS} ${MAKEFILE}
		@echo "Creating $(LIBRARY) ..."
		$(AR) cru $(LIBRARY) $(OBJS)
		$(RANLIB) $(LIBRARY)
		@echo "done"

clean:;		@rm -f $(OBJS) core

distclean:;	@rm -f $(OBJS) core $(LIBRARY) $(MAKEFILE) 

depend:;        @makedepend -f $(MAKEFILE) $(SRCS)

update:		$(DEST)/$(LIBRARY) 

# DO NOT DELETE THIS LINE -- make depend depends on it.
