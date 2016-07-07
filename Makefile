SHELL = /bin/sh

# List of subdirectories where files need to be compiled.
SUBDIRS = src_cpp examples_cpp

default: check all

# Check to make sure that the site configuration file is in place
check:
	@if [ ! -e "config/config.site" ]; then \
    echo "Please create the configuration file config/config.site first" ; \
    echo "See the README file in the directory config" ; \
    exit 1;\
  fi

# Generic make something command in a subdirectory
all clean test: check
	 for i in $(SUBDIRS) ;\
	 do echo "Making $@ in $$i ..." ;\
	   cd $$i; $(MAKE) $@ || exit 1; cd .. ;\
	 done

# Doxygen documentation
docs:
	 doxygen UQTkCppDoxyConfig
