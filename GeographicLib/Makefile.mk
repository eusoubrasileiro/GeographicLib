# $Id: Makefile.mk 6723 2009-10-18 21:33:42Z karney $

MAKEFILE := $(lastword $(MAKEFILE_LIST))
MAKE := $(MAKE) -f $(MAKEFILE)
SUBDIRS = src tools doc
ALLDIRS = include $(SUBDIRS) maxima windows

all: src tools

$(SUBDIRS):
	$(MAKE) -C $@

tools: src
install: install-headers install-lib install-tools
clean: clean-src clean-tools clean-doc
install-headers:
	$(MAKE) -C include install
install-lib:
	$(MAKE) -C src install
install-tools: src
	$(MAKE) -C tools install
install-doc: doc
	$(MAKE) -C doc install
clean-src:
	$(MAKE) -C src clean
clean-tools:
	$(MAKE) -C tools clean
clean-doc:
	$(MAKE) -C doc clean

list:
	@for f in 00README.txt COPYING AUTHORS Makefile $(MAKEFILE); do \
	  echo $$f; \
	done
	@for d in $(ALLDIRS); do \
	  (echo $(MAKEFILE); $(MAKE) -s -C $$d list) | tr ' ' '\n' | \
	  while read f; do echo $$d/$$f; done; \
	done

DATEVERSION:=$(shell date +%Y-%m)

package:
	test -d distrib || mkdir distrib
	$(MAKE) -s list | while read f;do \
	  echo GeographicLib/$$f; \
	done | xargs tar Ccfz .. distrib/Geographic-$(DATEVERSION).tgz
	rm -rf distrib/GeographicLib
	tar Cxfpz distrib distrib/Geographic-$(DATEVERSION).tgz
	find distrib/GeographicLib -type f \
	  \! \( -name "*.pdf" -o -name "*.png" \) | xargs unix2dos -k
	cd distrib && zip -r Geographic-$(DATEVERSION).zip GeographicLib && \
	  rm -rf GeographicLib

.PHONY: all $(SUBDIRS) \
	install install-headers install-lib install-tools install-doc \
	clean clean-src clean-tools clean-doc list package