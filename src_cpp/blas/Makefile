SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site
TARGET = libuqtkblas.a
F77FILES = $(wildcard *.f)
# F77FILES = caxpy.f ccopy.f cdotc.f cdotu.f cgbmv.f cgemm.f cgemv.f cgerc.f cgeru.f \
#            cscal.f csscal.f cswap.f ctbsv.f ctrmm.f ctrsm.f ctrsv.f d1mach.f dasum.f \
#            daxpy.f dcopy.f ddot.f dgbmv.f dgemm.f dgemv.f dger.f dnrm2.f drot.f drotg.f \
#            dscal.f dswap.f dsymv.f dsyr2.f dsyr2k.f dsyrk.f dtbsv.f dtrmm.f \
#            dtrmv.f dtrsm.f dzasum.f dznrm2.f icamax.f idamax.f isamax.f izamax.f sasum.f \
#            saxpy.f scasum.f scnrm2.f scopy.f sdot.f sgbmv.f sgemm.f sgemv.f sger.f \
#            snrm2.f srot.f srotg.f sscal.f sswap.f ssymv.f ssyr2.f stbsv.f strmm.f strsm.f \
#            zaxpy.f zcopy.f zdotc.f zdotu.f zdscal.f zgbmv.f zgemm.f zgemv.f zgerc.f zgeru.f \
#            zscal.f zswap.f ztbsv.f ztrmm.f ztrsm.f ztrsv.f
CXXFILES = 
CCFILES = 
INCFILES = blas.h
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
