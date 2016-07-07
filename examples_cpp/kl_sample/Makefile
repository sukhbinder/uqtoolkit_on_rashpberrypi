SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site

TARGET1 = kl_sample.x
SRCS1   = $(TARGET1:.x=.cpp)
OBJS1   = $(SRCS1:.cpp=.o)
INCS1   = $(SRCS1:.cpp=.h)

TARGET2 = kl_sample2D.x
SRCS2   = $(TARGET2:.x=.cpp)
OBJS2   = $(SRCS2:.cpp=.o)
INCS2   = $(SRCS2:.cpp=.h)

TARGET3 = kl_sample2Du.x
SRCS3   = $(TARGET3:.x=.cpp)
OBJS3   = $(SRCS3:.cpp=.o)
INCS3   = 

LIBINCDIR = $(TOPDIR)/src_cpp/include
LIBBINDIR = $(TOPDIR)/src_cpp/lib
INCDIRS = -I. -I$(LIBINCDIR)

LIBS = -L$(LIBBINDIR) -lUQTk -luqtktools -ldsfmt -l$(LAPACK) -l$(SLATEC) -l$(BLAS) $(FCLIB)
LDEP = $(LIBBINDIR)/libUQTk.a $(LIBBINDIR)/libuqtktools.a \
       $(LIBBINDIR)/libdsfmt.a $(LIBBINDIR)/lib$(LAPACK).a \
       $(LIBBINDIR)/lib$(SLATEC).a $(LIBBINDIR)/lib$(BLAS).a

DEFS = -D__$(FTNNAME)
.SUFFIXES : .cpp

all: $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET1): $(OBJS1) $(INCS1)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET1) $(OBJS1) $(LIBS)

$(TARGET2): $(OBJS2) $(INCS2)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET2) $(OBJS2) $(LIBS)

$(TARGET3): $(OBJS3) $(INCS3)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET3) $(OBJS3) $(LIBS)

clean:
	rm -f $(OBJS1) $(TARGET1) $(OBJS2) $(TARGET2) $(OBJS3) $(TARGET3) *~ *.pyc

.f.o:
	$(F77) $(FFLAGS) $(DEFS) $(INCDIRS) -c $*.f

.c.o:
	$(CC) $(CCFLAGS) $(DEFS) $(INCDIRS) -c $*.c

.cpp.o:
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -c $*.cpp

test: $(TARGET)
	@echo "The test target has not been implemented for $(TARGET)"

pres: doc/kl_example.tex
	cd doc; make; cd ..
