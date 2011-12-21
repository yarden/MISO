
all: Pythonpackage

########################################################

RSRC=$(wildcard src/*.c) $(wildcard src/*.h) $(wildcard src/*.pmt)
RSRC2=$(patsubst src/%,splicing/src/%,$(RSRC))
RSRC3=$(wildcard splicing/R/*.R) $(wildcard splicing/man/*.Rd) $(wildcard splicing/src/*.c) splicing/src/Makevars.in

Rpackage: splicing_0.1.tar.gz

splicing_0.1.tar.gz: $(RSRC2) $(RSRC3) splicing/DESCRIPTION splicing/NAMESPACE
	cd splicing && autoconf && autoheader
	R CMD build splicing

splicing/src/%.c: src/%.c
	cp src/$(@F) $@

splicing/src/%.h: src/%.h
	cp src/$(@F) $@

splicing/src/%.pmt: src/%.pmt
	cp src/$(@F) $@

tests: Rtests Pythontests

Rtests:
	cd splicing/tests && echo "tools:::.runPackageTestsR()" | \
        R --no-save && echo

Pythontests:
	cp pysplicing/test.py /tmp && cd /tmp && python test.py

########################################################

PSRC = $(wildcard src/*.c)
PSRC2 = $(wildcard src/*.h) $(wildcard src/*.pmt)
PSRC3 = $(patsubst src/%,pysplicing/src/%,$(PSRC))
PSRC4 = $(patsubst src/%,pysplicing/include/%,$(PSRC2))
PSRC5 = $(wildcard pysplicing/pysplicing/*.py)
PSRC6 = $(wildcard pysplicing/src/*.c) $(wildcard pysplicing/include/*.h) \
	$(wildcard pysplicing/src/lapack/*.c) \
	$(wildcard pysplicing/src/f2c/*.c)

Pythonpackage: pysplicing-0.1.tar.gz

pysplicing-0.1.tar.gz: $(PSRC3) $(PSRC4) $(PSRC5) $(PSRC6) pysplicing/setup.py pysplicing/MANIFEST.in
	rm -f pysplicing/MANIFEST
	cd pysplicing && python setup.py sdist -d ..

pysplicing/src/%.c: src/%.c
	cp src/$(@F) $@

pysplicing/include/%.h: src/%.h
	cp src/$(@F) $@

pysplicing/include/%.pmt: src/%.pmt
	cp src/$(@F) $@

.PHONY: all Rpackage tests Rtests Pythonpackage Pythontests
