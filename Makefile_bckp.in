# local variables
MAKE	      = make
SHELL         = /bin/sh
64_BIT_OS     = x86_64 alpha
MACHINE       = $(shell uname -m)

ifneq (,$(findstring $(MACHINE),$(64_BIT_OS)))
  FLAGS+=-D_64_BIT_OS_
endif

ifeq ($(MACHINE),alpha)
  FLAGS+=-D_NO_SOCKLEN_T_
endif

SUBDIRS	      = cmlib grid atmos slave jsub master # jstat jdel sdel

all:
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) all -C $$subdir FLAGS="$(FLAGS)" ; done
	@echo "done"

clean:
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) clean -C $$subdir TYPE_DEF="$(TYPE_DEF)" ; done

distclean: 
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) distclean -C $$subdir TYPE_DEF="$(TYPE_DEF)" ; done
	@rm -rf config.log config.h config.status Makefile autom4te.cache
dep:
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) dep -C $$subdir TYPE_DEF="$(TYPE_DEF)" ; done
install:
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) install -C $$subdir TYPE_DEF="$(TYPE_DEF)" ; done

makedepend:
	@makedepend -f$(MAKEFILE) $(OBJS:.o=.cc) $(INCLUDE)

update:
	$(DEST)/$(PROGRAM) 

# DO NOT DELETE THIS LINE -- make depend depends on it.
