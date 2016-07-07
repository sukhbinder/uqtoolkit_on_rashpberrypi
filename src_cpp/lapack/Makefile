SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site
TARGET = libuqtklapack.a
F77FILES = $(wildcard *.f)
# F77FILES = dgbtf2.f dgbtrf.f dgbtrs.f dgeqr2.f dgttrf.f dgttrs.f dlabad.f dlacon.f \
#            dlacpy.f dladiv.f dlae2.f dlaebz.f dlaev2.f dlaexc.f dlagtf.f dlagtm.f \
#            dlagts.f dlahqr.f dlaln2.f dlamch.f dlange.f dlanhs.f dlanst.f dlansy.f \
#            dlanv2.f dlaptm.f dlapy2.f dlapy3.f dlaran.f dlarfb.f dlarf.f dlarfg.f \
#            dlarft.f dlarfx.f dlarnd.f dlarnv.f dlartg.f dlaruv.f dlascl.f dlaset.f \
#            dlasr.f dlasrt.f dlassq.f dlaswp.f dlasy2.f dlatrd.f dorg2l.f dorg2r.f \
#            dorgql.f dorgqr.f dorgtr.f dorm2l.f dorm2r.f dormql.f dormqr.f dormtr.f \
#            dpttrf.f dpttrs.f dstebz.f dstein.f dsteqr.f dsterf.f dsyevx.f dsytd2.f \
#            dsytrd.f dtrevc.f dtrexc.f dtrsen.f dtrsyl.f dzsum1.f icmax1.f ilaenv.f \
#            izmax1.f lsame.f lsamen.f xerbla.f xlaenv.f dpotrf.f dpotf2.f dpotri.f \
#            dtrtri.f dlauum.f dlauu2.f dtrti2.f dgetrf.f dgetri.f dgetf2.f dgesv.f \
#            dgetrs.f dtrsm.f dgeev.f dgebak.f dgebal.f dgehrd.f dhseqr.f dorghr.f \
#            zgesv.f zgetrf.f zgetf2.f zgetrs.f zlaswp.f zgetri.f ztrti2.f ztrtri.f \
#            dgelq2.f dgeqrf.f dgelqf.f disnan.f dtrtrs.f dgels.f dorml2.f dlaisnan.f \
#            dormlq.f iladlc.f ieeeck.f iladlr.f iparmq.f lev2pchol.f dlaic1.f dstev.f
 
CXXFILES = 
CCFILES = 
INCFILES = lapack.h
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
