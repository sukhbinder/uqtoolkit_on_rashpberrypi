SHELL = /bin/sh

# List of subdirectories where files need to be compiled.
SUBDIRS = pce_sens pce_resp pce_quad generate_quad pce_rv pce_eval gen_mi

default: all


# Generic make something command in a subdirectory
all clean test links:
	 for i in $(SUBDIRS) ;\
	 do echo "Making $@ in $$i ..." ;\
	   cd $$i; $(MAKE) $@ || exit 1; cd .. ;\
	 done
