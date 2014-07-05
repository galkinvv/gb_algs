all:
	exec make SHELL=$(SHELL) BUILDDIR=/data/data/com.pdaxrom.cctools/root/cctools/home/build -j 4 check
c:
	exec make SHELL=$(SHELL) BUILDDIR=/data/data/com.pdaxrom.cctools/root/cctools/home/build clean