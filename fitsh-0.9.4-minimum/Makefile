SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld
INSTALL=install

.PHONY: all install \
        libfits src doc \
        info html pdf ps dvi \
        install install-strip \
        install-info install-html install-dvi install-pdf install-ps \
        help2man man \
        uninstall \
        clean distclean mostlyclean maintainter-clean TAGS \
        deb deb-fitsh

DEB_VERSION = 0.9.4
DEB_ARCH = amd64
DPKG_DEB = dpkg-deb
GZIP	 = gzip

datarootdir = ${prefix}/share
prefix = /usr/local
exec_prefix = ${prefix}
mandir = ${datarootdir}/man
bindir = ${exec_prefix}/bin
libdir = ${exec_prefix}/lib
includedir = ${prefix}/include

PACKAGE_NAME    = fitsh
PACKAGE_TARNAME = fitsh

infodir = ${datarootdir}/info
docdir  = ${datarootdir}/doc/${PACKAGE_TARNAME}
htmldir = ${docdir}
dvidir  = ${docdir}
pdfdir  = ${docdir}
psdir   = ${docdir}

DLEXT = so

HELP2MAN = true
ROFF	 = groff
MAN	 = man
DVIPS	 = dvips

ARCH_VERSION = $(DEB_VERSION)
DIST = fitsh-$(ARCH_VERSION)

DEB=misc/deb

all: src

libfits: 
	$(MAKE) -C libfits

src: libfits
	$(MAKE) -C src

doc:
	$(MAKE) -C doc

clean:
	$(MAKE) -C libfits clean
	$(MAKE) -C src clean
	$(MAKE) -C doc clean
	rm -f $(DIST).tar.gz
	rm -f *.deb

distclean: clean
	rm -f -r autom4te.cache
	rm -f ./config.status ./config.log
	rm -f ./Makefile
	rm -f ./config.h
	rm -f ./libfits/Makefile
	rm -f ./src/Makefile 
	rm -f ./src/index/Makefile 
	rm -f ./src/io/Makefile 
	rm -f ./src/link/Makefile 
	rm -f ./src/math/Makefile 
	rm -f ./src/math/dft/Makefile 
	rm -f ./src/math/elliptic/Makefile 
	rm -f ./src/math/elliptic/Makefile 
	rm -f ./src/math/expint/Makefile 
	rm -f ./src/math/fit/Makefile 
	rm -f ./src/math/intersec/Makefile 
	rm -f ./src/math/spline/Makefile 
	rm -f ./src/parallel/Makefile 
	rm -f ./doc/Makefile
	rm -f ./doc/examples/Makefile

install: src 
	$(INSTALL) -d $(DESTDIR)${bindir} 
	$(INSTALL) -m 0755 ./src/fistar $(DESTDIR)${bindir}
	$(INSTALL) -m 0755 ./src/grcollect $(DESTDIR)${bindir}
	$(INSTALL) -d $(DESTDIR)${libdir}/fitsh
	$(INSTALL) -d $(DESTDIR)${mandir}/man1
	$(INSTALL) -m 0644 ./man/fistar.1.gz $(DESTDIR)${mandir}/man1
	$(INSTALL) -m 0644 ./man/grcollect.1.gz $(DESTDIR)${mandir}/man1

uninstall:
	rm -f ${bindir}/fistar
	rm -f ${bindir}/grcollect
	rm -f ${mandir}/man1/fistar.1.gz
	rm -f ${mandir}/man1/grcollect.1.gz

HELP2MAN_LIST=fistar grcollect

help2man: src
	mkdir -p ./doc/help2man
	cat ./doc/man.brief | while read b name; do $(HELP2MAN) -n "$$name" -N -h --long-help ./src/$$b > ./doc/help2man/$$b.1; $(GZIP) -f ./doc/help2man/$$b.1; done

man: help2man
	mkdir	-p $(MAN)
	for b in $(HELP2MAN_LIST); do cp -p ./doc/help2man/$$b.1.gz $(MAN); done

dist:
	rm      -f -r $(DIST)
	mkdir   $(DIST)
	mkdir   $(DIST)/doc
	cp      -p doc/Makefile.in $(DIST)/doc
	mkdir	$(DIST)/doc/examples
	cp	-p doc/examples/Makefile.in $(DIST)/doc/examples
	mkdir   $(DIST)/include
	mkdir   $(DIST)/include/astro
	cp	-p include/astro/*.[ch] $(DIST)/include/astro
	mkdir   $(DIST)/include/fits
	cp	-p include/fits/*.[ch] $(DIST)/include/fits
	mkdir   $(DIST)/include/random
	cp	-p include/random/*.[ch] $(DIST)/include/random
	mkdir   $(DIST)/libfits
	cp      -p libfits/*.[ch] libfits/Makefile.in $(DIST)/libfits
	mkdir   $(DIST)/src
	cp      -p src/*.[ch] src/Makefile.in $(DIST)/src
	mkdir   $(DIST)/src/index
	cp      -p src/index/*.[ch] src/index/Makefile.in $(DIST)/src/index
	mkdir   $(DIST)/src/io
	cp      -p src/io/*.[ch] src/io/Makefile.in $(DIST)/src/io
	mkdir   $(DIST)/src/link
	cp      -p src/link/*.[ch] src/link/Makefile.in $(DIST)/src/link
	mkdir   $(DIST)/src/math
	cp      -p src/math/*.[ch] src/math/Makefile.in $(DIST)/src/math
	mkdir   $(DIST)/src/math/dft
	cp      -p src/math/dft/*.[ch] src/math/dft/Makefile.in $(DIST)/src/math/dft
	mkdir   $(DIST)/src/math/elliptic
	cp      -p src/math/elliptic/*.[ch] src/math/elliptic/Makefile.in $(DIST)/src/math/elliptic
	mkdir   $(DIST)/src/math/expint
	cp      -p src/math/expint/*.[ch] src/math/expint/Makefile.in $(DIST)/src/math/expint
	mkdir   $(DIST)/src/math/fit
	cp      -p src/math/fit/*.[ch] src/math/fit/Makefile.in $(DIST)/src/math/fit
	mkdir   $(DIST)/src/math/intersec
	cp      -p src/math/intersec/*.[ch] src/math/intersec/Makefile.in $(DIST)/src/math/intersec
	mkdir   $(DIST)/src/math/spline
	cp      -p src/math/spline/*.[ch] src/math/spline/Makefile.in $(DIST)/src/math/spline
	mkdir   $(DIST)/man
	cp      -p man/*.1.gz $(DIST)/man
	mkdir	$(DIST)/misc
	mkdir	$(DIST)/misc/deb
	cp	-p misc/deb/fitsh.control misc/deb/fitsh.list $(DIST)/misc/deb
	cp      -p AUTHORS CHANGELOG COPYING INSTALL README VERSION $(DIST)
	cp      -p ./configure configure.ac configure.mingw32 configure.mingw64 ./prepare $(DIST)
	cp      -p Makefile.in config.h.in $(DIST)
	tar     czf $(DIST).tar.gz $(DIST)
	rm      -f -r $(DIST)

deb: deb-fitsh

deb-fitsh: src help2man
	rm      -f -r $(DEB)/fitsh
	mkdir   -p $(DEB)/fitsh $(DEB)/fitsh/DEBIAN
	mkdir   -p $(DEB)/fitsh $(DEB)/fitsh/usr/bin $(DEB)/fitsh/usr/share/man/man1 $(DEB)/fitsh/usr/share/info $(DEB)/fitsh/usr/share/doc/fitsh
	cat $(DEB)/fitsh.list | while read file target ; do \
		cp $$file $(DEB)/fitsh/$$target ; \
	done
	cp      README  $(DEB)/fitsh/usr/share/doc/fitsh/
	gzip    $(DEB)/fitsh/usr/share/doc/fitsh/README
	cp      AUTHORS  $(DEB)/fitsh/usr/share/doc/fitsh/
	gzip    $(DEB)/fitsh/usr/share/doc/fitsh/AUTHORS
	cat     $(DEB)/fitsh.control | \
	sed     -e "s/__ARCH__/$(DEB_ARCH)/" \
	        -e "s/__VERSION__/$(DEB_VERSION)/" \
	        -e "s/__INSTALLED_SIZE__/`du -k -s --exclude DEBIAN $(DEB)/fitsh | awk '{ print $$1; }'`/" > $(DEB)/fitsh/DEBIAN/control
	for f in `find $(DEB)/fitsh` ; do \
	        o=`echo $$f | sed 's/deb\/fitsh\///'` ; \
	        if [ -f "$$f" -a `dirname "$$f"` != $(DEB)/fitsh/DEBIAN ] ; then \
	                m=`md5sum $$f | awk '{ print $$1; }'` ; \
	                echo $$m $$o ; \
	        fi ; \
	done > $(DEB)/fitsh/DEBIAN/md5sums
	$(DPKG_DEB) --build $(DEB)/fitsh fitsh_$(DEB_VERSION)_$(DEB_ARCH).deb
	rm      -f -r $(DEB)/fitsh
