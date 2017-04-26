all: build/libspooles.a

build/libspooles.a:
	rm -rf build ; mkdir build ; cd build ; tar xvfz ../archive/spooles.2.2.tgz
	cd build ; patch -p0 < ../patches/patch-ETree+src+makeGlobalLib 
	cd build ; patch -p0 < ../patches/patch-A2+makefile
	cd build ; patch -p0 < ../patches/patch-IVL+src+makeGlobalLib 
	cd build ; patch -p0 < ../patches/patch-MPI+makefile
	cd build ; patch -p0 < ../patches/patch-MPI+src+makeGlobalLib
	cd build ; patch -p0 < ../patches/patch-MT+drivers+AllInOneMT.c
	cd build ; patch -p0 < ../patches/patch-MT+makefile
	cd build ; patch -p0 < ../patches/patch-Make.inc
	cd build ; patch -p0 < ../patches/patch-Tree+src+makeGlobalLib
	cd build ; patch -p0 < ../patches/patch-Utilities+src+makeGlobalLib
	cd build ; patch -p0 < ../patches/patch-makefile
	cd build ; patch -p0 < ../patches/patch-timings.h
	cd build ; ${MAKE} global -f makefile
	cd build ; cp spooles.a libspooles.a
	
clean:
	rm -rf build

dist:
	
distdir:
	
