SUBDIRS =

if BUILD_GMP
SUBDIRS += gmp
endif

if BUILD_SPOOLES
SUBDIRS += spooles
endif

all:
	for i in $(SUBDIRS); do \
		echo "====> $(MAKE) $@ in $$i" ;\
		$(MAKE) -C $$i || exit 1; \
	done

clean:
	for i in $(SUBDIRS); do \
		echo "====> $(MAKE) $@ in $$i" ;\
		$(MAKE) -C $$i clean || exit 1; \
	done

