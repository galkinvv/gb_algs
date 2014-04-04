#!/usr/bin/make -f

ifndef CXX11
	CXX11=g++-4.8 -std=c++11
endif
ifndef WITH_MPI
	WITH_MPI=0
endif
ifndef OPTIMIZE
	OPTIMIZE = g
endif

GCC_WARNINGS=-Wall -Wextra -Wuninitialized -W -Wparentheses -Wformat=2 -Wswitch-default -Wcast-align -Wpointer-arith -Wwrite-strings -Wstrict-aliasing=2
GCC_WARNINGS_OFF=-Wno-missing-field-initializers -Wno-format-nonliteral -Wno-unknown-pragmas -Wno-reorder
ALL_CXX_LANG_FLAGS=-DWITH_MPI=$(WITH_MPI) $(GCC_WARNINGS_OFF) $(GCC_WARNINGS)

CXXFLAGS = $(ALL_CXX_LANG_FLAGS) -O$(OPTIMIZE) -ffunction-sections -fdata-sections -g -march=i686 -mtune=i686 -MD -MP
#-ffunction-sections -fdata-sections
#CXXFLAGS = -g -pg -O3 -march=i686 -mtune=i686

LDFLAGS = -O$(OPTIMIZE) -g 
#-Wl,--gc-sections
#-Wl,--gc-sections -Wl,--print-gc-sections
#LDFLAGS = -O3 -g -pg

MAINTARGET = integrtest/runalgo

BUILDDIR = build
ifeq ($(OS),Windows_NT)
ifeq ($(WITH_MPI),1)
CXXFLAGS+=-I"$(ProgramFiles)\MPICH2\include"
LDFLAGS+=-L"$(ProgramFiles)\MPICH2\lib" -lmpi
endif
RM=del
BINEXT=.exe
else
ifeq ($(WITH_MPI),1)
CXX11:=MPICH_CCC=$(firstword $(CXX11)) mpiCC $(wordlist 2,99,$(CXX11) )
endif
BINEXT=
RM=rm -rvf
endif
MAINBIN=$(MAINTARGET)$(BINEXT)

LD=$(CXX11)

MAINLIB=nmf5
FULLLIBNAME=$(BUILDDIR)/lib$(MAINLIB).a

LIBSOURCES = $(wildcard *.cpp)
ifeq ($(WITH_MPI),1)
	LIBSOURCES += $(wildcard mpi/*.cpp)
endif


ALLSOURCES=$(LIBSOURCES) $(MAINTARGET).cpp
LIBOBJECTS = $(LIBSOURCES:%.cpp=build/%.o)

all: $(BUILDDIR)/$(MAINBIN) $(BUILDDIR)/run-gt$(BINEXT)
$(BUILDDIR)/$(MAINBIN): $(BUILDDIR)/$(MAINTARGET).o $(FULLLIBNAME)
	$(LD) $(BUILDDIR)/$(MAINTARGET).o -L $(BUILDDIR) -l $(MAINLIB) -o $@ $(LDFLAGS)

TEST_SOURCES=$(wildcard unittest/*.cpp)
TEST_SOURCES+=$(wildcard unittest/mock/*.cpp)
TEST_OBJECTS = $(TEST_SOURCES:unittest/%.cpp=build/unittest/%.o)

$(BUILDDIR)/run-gt$(BINEXT): $(TEST_OBJECTS) $(FULLLIBNAME)
	$(LD) -pthread $^ -L $(BUILDDIR) -l $(MAINLIB) -o $@ $(LDFLAGS)

$(FULLLIBNAME): $(LIBOBJECTS)
	ar cr $@ $^

$(BUILDDIR)/unittest/%.o: unittest/%.cpp 3rd/gtest/src/gtest-all.cc
	mkdir -p $(dir $@)
	$(CXX11) $(CXXFLAGS) -I . -I 3rd/gtest -I 3rd/gtest/include -c $< -o $@

parse.tab:
	bison parse.ypp

$(BUILDDIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX11) $(CXXFLAGS) -I . -c $< -o $@

clean:
	$(RM) $(BUILDDIR)

quickcompile:
	$(CXX11) $(ALL_CXX_LANG_FLAGS) -I . -I 3rd/gtest -I 3rd/gtest/include -S -x c++ $(QUICK_SOURCE) -o /dev/null

3rd/gtest/src/gtest-all.cc:
	rm -rf 3rd/gtest/ /tmp/gtest.zip /tmp/gtest_version
	wget http://googletest.googlecode.com/files/gtest-1.7.0.zip -O /tmp/gtest.zip
	unzip /tmp/gtest.zip -d /tmp/gtest_version
	mkdir -p 3rd/
	mv /tmp/gtest_version/gtest* 3rd/gtest

3rd/gmp/include/gmp.h:
	rm -rf 3rd/gmp/ /tmp/gmp_build
	mkdir /tmp/gmp_build
	wget https://gmplib.org/download/gmp/gmp-6.0.0a.tar.xz -O /tmp/gmp_build/gmp-6.0.0a.tar.xz
	tar xf /tmp/gmp_build/gmp-6.0.0a.tar.xz -C /tmp/gmp_build
	mkdir -p 3rd/
	

check: $(BUILDDIR)/run-gt$(BINEXT)
	$(BUILDDIR)/run-gt$(BINEXT)

-include $(shell find $(BUILDDIR) -name '*.d')
