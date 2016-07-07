SHELL = /bin/sh
TOPDIR=../..
# Many macros are defined in the below included site-specific configuration file
include $(TOPDIR)/config/config.site

TARGET   = libuqtkmcmc.a
F77FILES = 
CXXFILES = uqtkmcmc.cpp
CCFILES  = 
INCFILES = uqtkmcmc.h
SRCS = $(CCFILES) $(CXXFILES) $(F77FILES)
OBJS = $(CCFILES:.c=.o) $(CXXFILES:.cpp=.o) $(F77FILES:.f=.o)

#TARGET1   = libtmcmc.a
#F77FILES1 = 
#CXXFILES1 = tmcmc.cpp
#CCFILES1  = 
#INCFILES1 = tmcmc.h
#SRCS1 = $(CCFILES1) $(CXXFILES1) $(F77FILES1)
#OBJS1 = $(CCFILES1:.c=.o) $(CXXFILES1:.cpp=.o) $(F77FILES1:.f=.o)

LIBINCDIR = ../include
LIBBINDIR = ../lib
INCDIRS = -I . -I $(LIBINCDIR)
DEFS = -D__$(FTNNAME) $(DFLAGS)
.SUFFIXES : .cpp


all: $(TARGET)
#all: $(TARGET) $(TARGET1)

$(TARGET): $(OBJS) links
	$(AR) $(TARGET) $(OBJS)
	$(RANLIB) $(TARGET)

#$(TARGET1): $(OBJS1) links
#	$(AR) $(TARGET1) $(OBJS1)
#	$(RANLIB) $(TARGET1)

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
#	for i in $(TARGET1) ;\
#	  do \rm -f $(LIBBINDIR)/$${i} ;\
#		ln -s $(PWD)/$${i} $(LIBBINDIR)/$${i} ;\
#	done
#	for i in $(INCFILES1) ;\
#	  do \rm -f $(LIBINCDIR)/$${i} ;\
#		ln -s $(PWD)/$${i} $(LIBINCDIR)/$${i} ;\
#	done

clean:
	rm -f $(OBJS) $(TARGET)
#clean:
#	rm -f $(OBJS) $(TARGET) $(OBJS1) $(TARGET1) 

.f.o:
	$(F77) $(FFLAGS) $(DEFS) $(INCDIRS) -c $*.f

.c.o:
	$(CC) $(CCFLAGS) $(DEFS) $(INCDIRS) -c $*.c

.cpp.o:
	$(CXX) $(CXXFLAGS) $(DEFS) $(INCDIRS) -c $*.cpp

test: $(TARGET)
	@echo "The test target has not been implemented for $(TARGET)"
