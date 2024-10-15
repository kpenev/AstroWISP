SHELL=/bin/sh

CC=gcc
AR=ar
LD=ld
RANLIB=ranlib

CFLAGS=-Wall -pedantic -ansi -O3 -Wno-strict-aliasing -fPIC -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_FITSH_SOURCE -Wno-overlength-strings -Wno-long-long -Wno-unused-result -I$(INC)

DEPTH=..

INC=$(DEPTH)/include

CFLARGEFILE=

DLEXT=so
DLSWC=-shared
DLDYN=-rdynamic
DLLIB=-ldl

TARGETS = fistar grcollect
TARGETS += 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

DEP_FIHEADER=\
	longhelp.o \
	str.o \
	ui.o

MOD_FIHEADER=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o

LIB_FIHEADER=$(DEPTH)/libfits/libfits.a -lm

DEP_FIINFO=\
	fitsmask.o \
	statistics.o \
	ui.o \
	longhelp.o \
	tensor.o

MOD_FIINFO=\
	math/fit/lmfit.o \
	math/poly.o \
	math/polyfit.o \
	math/splinefit.o \
	math/spline/biquad.o \
	math/spline/bicubic.o \
	math/spline/spline.o  \
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	link/link.a

LIB_FIINFO=$(DEPTH)/libfits/libfits.a -lm

DEP_FIPHOT=\
	fitsmask.o \
	statistics.o \
	magnitude.o \
	common.o \
	longhelp.o \
	ui.o \
	tensor.o \
	str.o \

MOD_FIPHOT=\
	math/intersec/intersec.o \
	math/intersec/intersec-cri.o \
	math/fit/lmfit.o \
	math/poly.o \
	math/spline/biquad.o \
	math/spline/biquad-isc.o \
	math/spline/spline.o \
	math/polygon.o \
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o

LIB_FIPHOT=$(DEPTH)/libfits/libfits.a -lm

DEP_FIRANDOM=\
	longhelp.o \
	ui.o \
	tensor.o \
	psf-io.o \
	star-draw.o \
	magnitude.o \
	common.o \
	fitsmask.o \

MOD_FIRANDOM=\
	math/fit/lmfit.o \
	math/poly.o \
	math/expint/expint.o \
	math/spline/biquad.o \
	math/spline/biquad-isc.o \
	math/intersec/intersec-cri.o \
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o

DEP_FISTAR=\
	fistar-io.o \
	fitsmask.o \
	statistics.o \
	tensor.o \
	stars.a \
	ui.o \
	longhelp.o \
	common.o \
	imgtrans.o \
	background.o \
	magnitude.o \
	psf.a

MOD_FISTAR=\
	math/fit/lmfit.o \
	math/poly.o \
	math/polyfit.o \
	math/spline/biquad.o \
	math/spline/biquad-isc.o \
	math/expint/expint.o \
	math/intersec/intersec-cri.o \
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	index/sort.o \
	index/multiindex.o \
	link/floodfill.o \
	link/link.a \
	math/delaunay.o

LIB_FISTAR=$(DEPTH)/libfits/libfits.a -lm

DEP_GRCOLLECT=\
	cache.o \
	statistics.o \
	longhelp.o \
	ui.o \
	str.o

MOD_GRCOLLECT=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o

LIB_GRCOLLECT=-lm

DEP_GRMATCH=\
	longhelp.o \
	ui.o \
	transform.o

MOD_GRMATCH=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	math/fit/lmfit.o \
	math/poly.o \
	math/polyfit.o \
	math/convexhull.o \
	math/match.a

LIB_GRMATCH=-lm

DEP_GRTRANS=\
	statistics.o \
	tensor.o \
	ui.o \
	longhelp.o \
	transform.o

MOD_GRTRANS=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	math/fit/lmfit.o \
	math/poly.o \
	math/polyfit.o

LIB_GRTRANS=-lm

DEP_GROPT=\
	longhelp.o \
	ui.o 

MOD_GROPT=\
	io/tokenize.o \
	io/scanarg.o \
	math/polygon.o \
	math/matrixvector.o 

LIB_GROPT=-lm

DEP_LFIT=\
	str.o \
	longhelp.o \
	ui.o \

MOD_LFIT=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	math/fit/lmfit.o \
	math/dft/pbfft.o \
	math/fit/downhill.o \
	math/elliptic/elliptic.o \
	math/elliptic/ntiq.o \
	math/spline/spline.o

#   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #

DEP_FIC_MPSTACK=\
	fitsmask.o \
	statistics.o \
	transform.o \
	common.o \
	ui.o \
	longhelp.o \
	tensor.o

MOD_FIC_MPSTACK=\
	io/iof.o \
	io/scanarg.o \
	io/tokenize.o \
	math/spline/bicubic.o \
	math/spline/spline.o \
	math/poly.o

LIB_FIC_MPSTACK=$(DEPTH)/libfits/libfits.a -lm

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

DIRS=index io link math

all: dirs
	$(MAKE) targets

targets: $(TARGETS)

dirs: $(DIRS)

.PHONY: all clean dirs $(DIRS) targets

###############################################################################

index:
	$(MAKE) -C index

io:
	$(MAKE) -C io

link:
	$(MAKE) -C link

math:
	$(MAKE) -C math

###############################################################################

background.o: background.c background.h
	$(CC) $(CFLAGS) -c background.c

cache.o: cache.c cache.h
	$(CC) $(CFLAGS) -c cache.c

imgtrans.o: imgtrans.c imgtrans.h
	$(CC) $(CFLAGS) -c imgtrans.c

magnitude.o: magnitude.c magnitude.h
	$(CC) $(CFLAGS) -c magnitude.c

statistics.o: statistics.c statistics.h
	$(CC) $(CFLAGS) -c statistics.c

str.o: str.c str.h
	$(CC) $(CFLAGS) -c str.c

tensor.o: tensor.c tensor.h
	$(CC) $(CFLAGS) -c tensor.c

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

fitsmask.o: fitsmask.c fitsmask.h
	$(CC) $(CFLAGS) -c fitsmask.c

common.o: common.c common.h fitsh.h
	$(CC) $(CFLAGS) -c common.c

psf-base.o: psf-base.c psf-base.h psf.h
	$(CC) $(CFLAGS) -c psf-base.c

psf-io.o: psf-io.c psf-io.h psf.h
	$(CC) $(CFLAGS) -c psf-io.c

psf-determine.o: psf-determine.c psf-determine.h psf.h stars.h
	$(CC) $(CFLAGS) -c psf-determine.c

ui.o: ui.c fitsh.h
	$(CC) $(CFLAGS) -c ui.c

star-base.o: star-base.c stars.h
	$(CC) $(CFLAGS) -c star-base.c

star-draw.o: star-draw.c stars.h
	$(CC) $(CFLAGS) -c star-draw.c

star-model.o: star-model.c star-model.h stars.h
	$(CC) $(CFLAGS) -c star-model.c

star-psf.o: star-psf.c stars.h
	$(CC) $(CFLAGS) -c star-psf.c

star-cand-pp.o: star-cand-pp.c stars.h
	$(CC) $(CFLAGS) -c star-cand-pp.c

star-cand-biq.o: star-cand-biq.c stars.h
	$(CC) $(CFLAGS) -c star-cand-biq.c

star-cand-trb.o: star-cand-trb.c stars.h
	$(CC) $(CFLAGS) -c star-cand-trb.c

star-cand-lnk.o: star-cand-lnk.c stars.h
	$(CC) $(CFLAGS) -c star-cand-lnk.c

MOD_STARS=star-base.o star-model.o star-psf.o star-draw.o \
          star-cand-pp.o star-cand-biq.o star-cand-trb.o star-cand-lnk.o 

MOD_PSF=psf-io.o psf-base.o psf-determine.o

stars.a: $(MOD_STARS)
	$(AR) src stars.a $(MOD_STARS)
	$(RANLIB) stars.a

psf.a: $(MOD_PSF)
	$(AR) src psf.a $(MOD_PSF)
	$(RANLIB) psf.a

longhelp.o: longhelp.c longhelp.h
	$(CC) $(CFLAGS) -c longhelp.c

###############################################################################

fistar-io.o: fistar-io.c fistar.h
	$(CC) $(CFLAGS) -c fistar-io.c

fistar: fistar.c fitsh.h $(DEP_FISTAR) $(MOD_FISTAR) fistar.h
	$(CC) $(CFLAGS) -o fistar fistar.c $(DEP_FISTAR) $(MOD_FISTAR) $(LIB_FISTAR)

grcollect: grcollect.c fitsh.h $(DEP_GRCOLLECT) $(MOD_GRCOLLECT)
	$(CC) $(CFLAGS) -o grcollect grcollect.c $(DEP_GRCOLLECT) $(MOD_GRCOLLECT) $(LIB_GRCOLLECT)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

fic_mpstack: fic_mpstack.c fitsh.h $(DEP_FIC_MPSTACK)
	$(CC) $(CFLAGS) -o fic_mpstack fic_mpstack.c $(DEP_FIC_MPSTACK) $(MOD_FIC_MPSTACK) $(LIB_FIC_MPSTACK)

###############################################################################

clean:
	$(MAKE) -C index clean
	$(MAKE) -C io clean
	$(MAKE) -C link clean
	$(MAKE) -C math clean
	rm -f *.o *.a *.exe $(TARGETS)

