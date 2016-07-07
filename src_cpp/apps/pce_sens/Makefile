SHELL = /bin/sh
TOPDIR=../../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site
TARGET = pce_sens
SRCS = pce_sens.cpp 
OBJS = $(SRCS:.cpp=.o)
INCS = 
LIBINCDIR = $(TOPDIR)/src_cpp/include
LIBDIR    = $(TOPDIR)/src_cpp/lib
LIBBINDIR = $(TOPDIR)/src_cpp/bin
INCDIRS = -I . -I $(LIBINCDIR)

LIBS = -L$(LIBDIR) -lUQTk -lquad -luqtktools -lcvode-2.6.0 -ldsfmt -l$(LAPACK) -l$(SLATEC) -l$(BLAS) $(FCLIB) 
LDEP = $(LIBDIR)/libUQTk.a $(LIBDIR)/libquad.a $(LIBDIR)/libuqtktools.a $(LIBDIR)/libcvode-2.6.0.a \
       $(LIBDIR)/libdsfmt.a $(LIBDIR)/lib$(LAPACK).a $(LIBDIR)/lib$(SLATEC).a $(LIBDIR)/lib$(BLAS).a
DEFS = -D__$(FTNNAME)
.SUFFIXES : .cpp

all: $(TARGET)

$(TARGET): $(OBJS) $(INCS) $(LDEP) links
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET) $(OBJS) $(LIBS)

# Put symbolic links to the executable files in the designated binary directory
links: 
	if ! [ -d $(LIBBINDIR) ] ; then \rm -f $(LIBBINDIR); mkdir $(LIBBINDIR); fi
	for i in $(TARGET) ;\
	  do \rm -f $(LIBBINDIR)/$${i} ;\
		ln -s $(PWD)/$${i} $(LIBBINDIR)/$${i} ;\
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
