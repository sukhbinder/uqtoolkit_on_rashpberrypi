SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site
TARGET = libuqtkslatec.a
F77FILES = fdump.f  i1mach.f  j4save.f  xercnt.f  xerhlt.f  xermsg.f  xerprn.f  xersve.f  xgetua.f \
           dasum.f dgbfa.f dgbsl.f dgeco.f dgefa.f dgefs.f dgesl.f dgmres.f dhels.f dheqr.f dorth.f dpigmr.f \
           drlcal.f dsort.f dxlcal.f isdgmr.f isort.f 
CXXFILES = 
CCFILES = 
INCFILES = slatec.h
SRCS = $(CCFILES) $(CXXFILES) $(F77FILES)
OBJS = $(CCFILES:.c=.o) $(CXXFILES:.cpp=.o) $(F77FILES:.f=.o)
LIBINCDIR = ../include
LIBBINDIR = ../lib
INCDIRS = -I . -I $(LIBINCDIR)
DEFS = -D__$(FTNNAME)
.SUFFIXES : .cpp


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

.cpp.o:
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -c $*.cpp

test: $(TARGET)
	@echo "The test target has not been implemented for $(TARGET)"
