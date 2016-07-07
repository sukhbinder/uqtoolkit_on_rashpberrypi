SHELL = /bin/sh

# List of subdirectories where files need to be compiled.
SUBDIRS = line_infer surf_rxn ops pce_bcs test_quad

default: all


# Generic make something command in a subdirectory
all clean test:
	 for i in $(SUBDIRS) ;\
	 do echo "Making $@ in $$i ..." ;\
	   cd $$i; $(MAKE) $@ || exit 1; cd .. ;\
	 done
