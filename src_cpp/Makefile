SHELL = /bin/sh
LIBDIR=lib
BINDIR=bin
INCDIR=include

# List of subdirectories where files need to be compiled.
SUBDIRS = slatec lapack blas lbfgs utils uqtktools array dsfmt quad cvode-2.6.0 uqtk bcs uqtkmcmc xmlutils apps

default: all

# make all (depends on links so that all .h file links are in place)
all: links
	 for i in $(SUBDIRS) ;\
	 do echo "Making $@ in $$i ..." ;\
	   cd $$i; $(MAKE) $@ || exit 1; cd .. ;\
	 done

# Generic make something command in a subdirectory
links test:
	 for i in $(SUBDIRS) ;\
	 do echo "Making $@ in $$i ..." ;\
	   cd $$i; $(MAKE) $@ || exit 1; cd .. ;\
	 done

# make clean also cleans out lib, bin, and include directories
clean:
	 for i in $(SUBDIRS) ;\
	 do echo "Making $@ in $$i ..." ;\
	   cd $$i; $(MAKE) $@; cd .. ;\
	 done
	 \rm -rf $(LIBDIR)
	 \rm -rf $(BINDIR)
	 \rm -rf $(INCDIR)
