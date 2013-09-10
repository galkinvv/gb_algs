#!/usr/bin/make -f

WITH_MPI=0
OPTIMIZE = 3
CXXFLAGS = -O$(OPTIMIZE) -g -march=i686 -mtune=i686 -std=c++0x -DWITH_MPI=$(WITH_MPI)
#-ffunction-sections -fdata-sections
#CXXFLAGS = -g -pg -O3 -march=i686 -mtune=i686

LDFLAGS = -O$(OPTIMIZE) -g
#-Wl,--gc-sections -Wl,--print-gc-sections
#LDFLAGS = -O3 -g -pg

MAINTARGET = testnmf5

BUILDDIR = build
ifeq ($(OS),Windows_NT)
ifeq ($(WITH_MPI),1)
CXXFLAGS+=-I"$(ProgramFiles)\MPICH2\include"
LDFLAGS+=-L"$(ProgramFiles)\MPICH2\lib" -lmpi
endif
CXX=g++
RM=del
BINEXT=.exe
else
ifeq ($(WITH_MPI),1)
CXX=mpiCC
else
CXX=g++
endif
BINEXT=
RM=rm -f
endif
MAINBIN=$(MAINTARGET)$(BINEXT)

LD=$(CXX)

MAINLIB=nmf5
FULLLIBNAME=$(BUILDDIR)/lib$(MAINLIB).a

LIBSOURCES = \
	cmatrix.cpp\
	cmodular.cpp\
	cmonomial.cpp\
	commonpolyops.cpp\
	cpolynomial.cpp\
	f4main.cpp\
	gbimpl.cpp\
	globalf4.cpp\
	libf4mpi.cpp\
	monomialmap.cpp\
	outputroutines.cpp\
	parse.tab.cpp\
	f5c_plain.cpp\
	f5_plain.cpp\
	reducebyset.cpp


ifeq ($(WITH_MPI),1)
	LIBSOURCES += mpimatrix.cpp
endif


ALLSOURCES=$(LIBSOURCES) $(MAINTARGET).cpp
LIBOBJECTS = $(LIBSOURCES:%.cpp=build/%.o)
ALLOBJECTS = $(ALLSOURCES:%.cpp=build/%.o)

$(MAINBIN): $(BUILDDIR)/$(MAINTARGET).o $(FULLLIBNAME)
	$(LD) $(BUILDDIR)/$(MAINTARGET).o -L $(BUILDDIR) -l $(MAINLIB) -o $@ $(LDFLAGS)

TEST_SOURCES=tst/run-gt.cpp\
	tst/cmodular.cpp
TEST_OBJECTS = $(TEST_SOURCES:tst/%.cpp=build/tst/%.o)

run-gt$(BINEXT): $(TEST_OBJECTS) $(FULLLIBNAME)
	$(LD) -pthread $^ -L $(BUILDDIR) -l $(MAINLIB) -o $@ $(LDFLAGS)

$(FULLLIBNAME): $(LIBOBJECTS)
	ar cr $@ $^

$(BUILDDIR)/tst/run-gt.o: 3rd/gtest/src/gtest-all.cc

$(BUILDDIR)/tst/%.o: tst/%.cpp
	mkdir -p $(BUILDDIR)/tst
	$(CXX) $(CXXFLAGS)  -I . -I 3rd/gtest -I 3rd/gtest/include -c $< -o $@

$(BUILDDIR)/%.o: %.cpp
	mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(MAINTARGET) $(ALLOBJECTS) $(MAINBIN) $(FULLLIBNAME)

3rd/gtest/src/gtest-all.cc:
	rm -rf 3rd/gtest/ /tmp/gtest.zip /tmp/gtest_version
	wget http://googletest.googlecode.com/files/gtest-1.7.0-rc1.zip -O /tmp/gtest.zip
	unzip /tmp/gtest.zip -d /tmp/gtest_version
	mkdir -p 3rd/
	mv /tmp/gtest_version/gtest* 3rd/gtest