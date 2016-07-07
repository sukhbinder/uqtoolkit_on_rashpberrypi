SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site
TARGET = libcvode-2.6.0.a
F77FILES = 
CXXFILES = 
CCFILES = cvode.c cvode_band.c cvode_bandpre.c cvode_bbdpre.c cvode_dense.c \
          cvode_diag.c cvode_direct.c cvode_io.c cvode_lapack.c cvode_spbcgs.c \
          cvode_spgmr.c cvode_spils.c cvode_sptfqmr.c \
          sundials_band.c sundials_dense.c sundials_direct.c sundials_iterative.c \
          sundials_math.c sundials_nvector.c sundials_spbcgs.c sundials_spgmr.c \
          sundials_sptfqmr.c nvector_serial.c 
INCCVFILES = cvode_band.h cvode_bandpre_impl.h cvode_bbdpre_impl.h cvode_diag.h \
             cvode_direct.h cvode.h cvode_lapack.h cvode_spgmr.h cvode_spils_impl.h \
             cvode_bandpre.h cvode_bbdpre.h cvode_dense.h cvode_diag_impl.h \
             cvode_direct_impl.h cvode_impl.h cvode_spbcgs.h cvode_spils.h cvode_sptfqmr.h
INCNVFILES = nvector_parallel.h  nvector_serial.h
INCSDFILES = sundials_band.h sundials_dense.h sundials_fnvector.h sundials_lapack.h \
             sundials_nvector.h sundials_spgmr.h sundials_types.h sundials_config.h sundials_direct.h \
             sundials_iterative.h sundials_math.h sundials_spbcgs.h sundials_sptfqmr.h

#INCFILES = cvode.h cvode_band.h cvode_bandpre.h cvode_bandpre_impl.h cvode_bbdpre.h \
#           cvode_bbdpre_impl.h cvode_dense.h cvode_diag.h cvode_diag_impl.h         \
#           cvode_direct.h cvode_direct_impl.h cvode_impl.h cvode_lapack.h           \
#           cvode_spbcgs.h cvode_spgmr.h cvode_spils.h cvode_spils_impl.h cvode_sptfqmr.h \
#           sundials_band.h sundials_config.h sundials_dense.h sundials_direct.h          \
#           sundials_fnvector.h sundials_iterative.h sundials_lapack.h sundials_math.h    \
#           sundials_nvector.h sundials_spbcgs.h sundials_spgmr.h sundials_sptfqmr.h      \
#           sundials_types.h nvector_serial.h 
SRCS = $(CCFILES) $(CXXFILES) $(F77FILES)
OBJS = $(CCFILES:.c=.o) $(CXXFILES:.cpp=.o) $(F77FILES:.f=.o)
LIBINCDIR = ../include
LIBBINDIR = ../lib
INCDIRS = -I. -I$(LIBINCDIR)
DEFS = -D__$(FTNNAME)
.SUFFIXES : .cpp


all: $(TARGET)

$(TARGET): $(OBJS) links
	$(AR) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

$(OBJS): links

# Put symbolic links to the library and include files in the designated binary and include directories
# This is slightly different from other libraries as cvode expects to see headers in specific subdirectories
links: 
	if ! [ -d $(LIBBINDIR) ] ; then \rm -f $(LIBBINDIR); mkdir $(LIBBINDIR); fi
	for i in $(TARGET) ;\
	  do \rm -f $(LIBBINDIR)/$${i} ;\
		ln -s $(PWD)/$${i} $(LIBBINDIR)/$${i} ;\
	done
	if ! [ -d $(LIBINCDIR) ] ; then \rm -f $(LIBINCDIR); mkdir $(LIBINCDIR); fi
	if ! [ -d $(LIBINCDIR)/cvode ] ; then mkdir $(LIBINCDIR)/cvode; fi
	if ! [ -d $(LIBINCDIR)/nvector ] ; then mkdir $(LIBINCDIR)/nvector; fi
	if ! [ -d $(LIBINCDIR)/sundials ] ; then mkdir $(LIBINCDIR)/sundials; fi
	for i in $(INCCVFILES) ;\
          do \rm -f $(LIBINCDIR)/cvode/$${i} ;\
                ln -s $(PWD)/$${i} $(LIBINCDIR)/cvode/$${i} ;\
        done
	for i in $(INCNVFILES) ;\
          do \rm -f $(LIBINCDIR)/nvector/$${i} ;\
                ln -s $(PWD)/$${i} $(LIBINCDIR)/nvector/$${i} ;\
        done
	for i in $(INCSDFILES) ;\
          do \rm -f $(LIBINCDIR)/sundials/$${i} ;\
                ln -s $(PWD)/$${i} $(LIBINCDIR)/sundials/$${i} ;\
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
