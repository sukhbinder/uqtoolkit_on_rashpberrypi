SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site
TARGET = test_quad.x
SRCS   = test_quad.cpp
OBJS   = $(SRCS:.cpp=.o)
LIBINCDIR = $(TOPDIR)/src_cpp/include
LIBBINDIR = $(TOPDIR)/src_cpp/lib
INCDIRS = -I. -I$(LIBINCDIR)
LIBS = -L$(LIBBINDIR) -luqtktools -l$(LAPACK) -l$(SLATEC) -l$(BLAS) $(FCLIB) 
LDEP = $(LIBBINDIR)/libuqtktools.a $(LIBBINDIR)/lib$(LAPACK).a $(LIBBINDIR)/lib$(SLATEC).a $(LIBBINDIR)/lib$(BLAS).a
DEFS = -D__$(FTNNAME) $(DFLAGS)
.SUFFIXES : .cpp

all: $(TARGET)

$(TARGET): $(OBJS) $(LDEP)
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -o $(TARGET) $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)

.f.o:
	$(F77) $(FFLAGS) $(DEFS) $(INCDIRS) -c $*.f

.c.o:
	$(CC) $(CCFLAGS) $(DEFS) $(INCDIRS) -c $*.c

.cpp.o:
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -c $*.cpp

test: $(TARGET)
	@./$(TARGET); ./compare.x
