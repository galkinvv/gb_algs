#!/usr/bin/make -f

ifndef CXX11
	CXX11=g++-4.8
endif
ifndef CC_FORCXX
	CC_FORCXX=gcc-4.8
endif
ifndef WITH_MPI
	WITH_MPI=0
endif
ifndef OPTIMIZE
	OPTIMIZE = 0
endif
ifndef BUILDDIR
	BUILDDIR = build
endif

GCC_WARNINGS=-Wall -Wextra -Wuninitialized -W -Wparentheses -Wformat=2 -Wswitch-default -Wcast-align -Wpointer-arith -Wwrite-strings -Wstrict-aliasing=2
GCC_WARNINGS_OFF=-Wno-missing-field-initializers -Wno-format-nonliteral -Wno-unknown-pragmas -Wno-reorder
ALL_CXX_LANG_FLAGS=-DWITH_MPI=$(WITH_MPI) $(GCC_WARNINGS_OFF) $(GCC_WARNINGS) -std=c++11

CXXFLAGS = $(ALL_CXX_LANG_FLAGS) -O$(OPTIMIZE) -g -march=native -mtune=native -MD -MP -ffunction-sections -fdata-sections
#CXXFLAGS = -g -pg -O3 -march=native -mtune=native

LDFLAGS = -O$(OPTIMIZE) -g 
#-Wl,--gc-sections
#-Wl,--gc-sections -Wl,--print-gc-sections
#LDFLAGS = -O3 -g -pg

MAINTARGET = runalgo

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

LD=$(CXX11) -std=c++11

MAINLIB=nmf5
OBJDIR=$(BUILDDIR)/obj
TMPDIR=$(BUILDDIR)/tmp
FULLLIBNAME=$(OBJDIR)/lib$(MAINLIB).a

LIBSOURCES = $(wildcard *.cpp)
ifeq ($(WITH_MPI),1)
	LIBSOURCES += $(wildcard mpi/*.cpp)
endif


LIBOBJECTS = $(LIBSOURCES:%.cpp=$(OBJDIR)/%.o)

all: $(BUILDDIR)/$(MAINBIN) $(BUILDDIR)/run-gt$(BINEXT)
$(BUILDDIR)/$(MAINBIN): $(OBJDIR)/testapps/$(MAINTARGET).o $(FULLLIBNAME)
	$(LD) $< -L $(OBJDIR) -l $(MAINLIB) -L 3rd/gmp/lib -l gmp -l gmpxx -o $@ $(LDFLAGS)

TEST_SOURCES=$(wildcard libtests/*.cpp)
TEST_SOURCES+=$(wildcard libtests/mock/*.cpp)
TEST_OBJECTS = $(TEST_SOURCES:libtests/%.cpp=$(OBJDIR)/libtests/%.o)

$(BUILDDIR)/run-gt$(BINEXT): $(TEST_OBJECTS) $(FULLLIBNAME)
	$(LD) -pthread $^ -L $(OBJDIR) -l $(MAINLIB) -L 3rd/gmp/lib -l gmp -l gmpxx -o $@ $(LDFLAGS)

$(FULLLIBNAME): $(LIBOBJECTS)
	ar cr $@ $^

$(OBJDIR)/libtests/%.o: libtests/%.cpp 3rd/gtest/src/gtest-all.cc 3rd/gmp/include/gmp.h
	mkdir -p $(dir $@)
	$(CXX11) $(CXXFLAGS) -I . -I 3rd/gtest -I 3rd/gtest/include -I 3rd/gmp/include/ -c $< -o $@

$(OBJDIR)/%.o: %.cpp 3rd/gmp/include/gmp.h
	mkdir -p $(dir $@)
	$(CXX11) $(CXXFLAGS) -I . -I 3rd/gmp/include/ -c $< -o $@

parse.tab:
	bison parse.ypp

clean:
	$(RM) $(BUILDDIR)
FD=$(shell pwd)
QUICKCOMPILE_OPTIONS = $(ALL_CXX_LANG_FLAGS) -I $(FD) -I $(FD)/3rd/gtest -I $(FD)/3rd/gtest/include -I $(FD)/3rd/gmp/include/ -S -o /dev/null

quickcompile_options: 
	echo $(QUICKCOMPILE_OPTIONS) 

quickcompile: 3rd/gtest/src/gtest-all.cc 3rd/gmp/include/gmp.h
	$(CXX11) -x c++ $(QUICKCOMPILE_OPTIONS) $(QUICK_SOURCE)

3rd/gtest/src/gtest-all.cc:
	rm -rf 3rd/gtest/ $(TMPDIR)/gtest.zip tmp/gtest_build
	mkdir -p $(TMPDIR)/gtest_build
	wget http://googletest.googlecode.com/files/gtest-1.7.0.zip -O $(TMPDIR)/gtest.zip
	unzip $(TMPDIR)/gtest.zip -d $(TMPDIR)/gtest_build
	mkdir -p 3rd/
	mv $(TMPDIR)/gtest_build/gtest* 3rd/gtest

3rd/gmp/include/gmp.h:
	rm -rf 3rd/gmp/ $(TMPDIR)/gmp_build
	mkdir -p $(TMPDIR)/gmp_build
	wget https://gmplib.org/download/gmp/gmp-6.0.0a.tar.xz -O $(TMPDIR)/gmp.tar.xz
	tar xf $(TMPDIR)/gmp.tar.xz -C $(TMPDIR)/gmp_build
	mkdir -p 3rd/
	cd $(TMPDIR)/gmp_build/gmp-* && CC=$(CC_FORCXX) CXX=$(firstword $(CXX11)) CXXFLAGS=$(wordlist 2,999,$(CXX11)) ./configure --prefix=$(TMPDIR)/gmp_install --enable-cxx=yes && make -j 4 && make install
	cp -r $(TMPDIR)/gmp_install 3rd/gmp

check: $(BUILDDIR)/run-gt$(BINEXT)
	$^

-include $(shell find $(BUILDDIR) -name '*.d')
