SPOOLESVERSION=2.2

EXTRACT_STAMP=.spooles_extract_done
PATCH_STAMP=.spooles_patch_done
CONFIGURE_STAMP=.spooles_configure_done
BUILD_STAMP=.spooles_build_done
INSTALL_STAMP=.spooles_install_done
INSTALL_INTERNAL_STAMP=.spooles_install_internal_done
CHECK_STAMP=.spooles_check_done

all: extract config build check install_internal

work/$(EXTRACT_STAMP):
	rm -rf work ; mkdir -p work/internal
	cd $(abs_builddir)/work/internal; tar xvfz ../../download/spooles.$(SPOOLESVERSION).tgz
	touch $(abs_builddir)/work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-A2+makefile
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-ETree+src+makeGlobalLib
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-IVL+src+makeGlobalLib
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-MPI+makefile
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-MPI+src+makeGlobalLib
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-MT+drivers+AllInOneMT.c
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-MT+makefile
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-Make.inc
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-Tree+src+makeGlobalLib
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-Utilities+src+makeGlobalLib
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-makefile
	cd $(abs_builddir)/work/internal; patch -p0 < ../../patches/patch-timings.h
	touch $(abs_builddir)/work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
	touch $(abs_builddir)/work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd $(abs_builddir)/work/internal/ ; ${MAKE} global -f makefile
#	cd $(abs_builddir)/work/internal/MT/src ; ${MAKE} global -f makeGlobalLib
	touch $(abs_builddir)/work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
#	cd $(abs_builddir)/work/internal/ ; ${MAKE} check
	touch $(abs_builddir)/work/$(CHECK_STAMP)

work/$(INSTALL_INTERNAL_STAMP): work/$(CHECK_STAMP)
	rm -rf ../i/SPOOLES/lib
	mkdir -p ../i/SPOOLES/lib
	cd $(abs_builddir)/work/internal/    ; cp spooles.a ../../../i/SPOOLES/lib/libspooles.a
	touch $(abs_builddir)/work/$(INSTALL_INTERNAL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install_internal: work/$(INSTALL_INTERNAL_STAMP)

distdir:

clean:
	rm -rf work
distclean:
	rm -rf work
