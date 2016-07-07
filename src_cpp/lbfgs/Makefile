SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site
TARGET = liblbfgs.a
F77FILES = lbfgs_routines.f
CCFILES  = lbfgsDR.c
INCFILES = lbfgs_routines.h

SRCS = $(CCFILES) $(F77FILES)
OBJS = $(CCFILES:.c=.o) $(F77FILES:.f=.o)
LIBINCDIR = ../include
LIBBINDIR = ../lib
INCDIRS = -I. -I$(LIBINCDIR)
DEFS = -D__$(FTNNAME) 
.SUFFIXES : .c .f


all: $(TARGET)

$(TARGET): $(OBJS) links
	$(AR) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

# Put symbolic links to the library and include files in the designated binary and include directories
links: 
	if ! [ -d $(LIBBINDIR) ] ; then \rm -f $(LIBBINDIR); mkdir $(LIBBINDIR); fi
	for i in $(TARGET) ;\
	  do \rm -f $(LIBBINDIR)/$${i} ;\
		ln -s $(PWD)/$${i} $(LIBBINDIR)/$${i} ;\
	done
	if ! [ -d $(LIBINCDIR) ] ; then \rm -f $(LIBINCDIR); mkdir $(LIBINCDIR); fi
	for i in $(INCFILES) ;\
	  do \rm -f $(LIBINCDIR)/$${i} ;\
		ln -s $(PWD)/$${i} $(LIBINCDIR)/$${i} ;\
	done

clean:
	rm -f $(OBJS) $(TARGET) 

.f.o:
	$(F77) $(FFLAGS) $(DEFS) $(INCDIRS) -c $*.f

.c.o:
	$(CC) $(CCFLAGS) $(DEFS) $(INCDIRS) -c $*.c

test: $(TARGET)
	@echo "The test target has not been implemented for $(TARGET)"
