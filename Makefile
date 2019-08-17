ROOTDIR = $(CURDIR)

ifeq ($(OS), Windows_NT)
	UNAME="Windows"
else
	UNAME=$(shell uname)
endif

# set compiler defaults for OSX versus *nix
# let people override either
OS := $(shell uname)
ifeq ($(OS), Darwin)
ifndef CC
export CC = $(if $(shell which clang), clang, gcc)
endif
ifndef CXX
export CXX = $(if $(shell which clang++), clang++, g++)
endif
else
# linux defaults
ifndef CC
export CC = gcc
endif
ifndef CXX
export CXX = g++
endif
endif

export LDFLAGS = 
export CFLAGS = -std=c++11 -Wall -Wno-unknown-pragmas -I ./include 

ifneq ($(UNAME), Windows)
	CFLAGS += -fPIC
ifeq ($(UNAME), Linux)
	LDFLAGS += -lrt
	PSM_DYLIB = lib/libpsm.so
else
	PSM_DYLIB = lib/libpsm.dylib
endif
else
	PSM_DYLIB = lib/libpsm.dll
endif

# specify tensor path
.PHONY: clean all clean_all doxygen Pypack Pyinstall Rpack Rbuild Rcheck


all: lib/libpsm.a $(PSM_DYLIB)

dylib: $(PSM_DYLIB)

build/PSM.o: src/PSM.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CFLAGS) $< -o $@

build/api.o: src/api.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CFLAGS) $< -o $@

lib/libpsm.a: build/PSM.o build/api.o
	@mkdir -p $(@D)
	ar crv $@ $^

lib/libpsm.dll lib/libpsm.so lib/libpsm.dylib: build/PSM.o build/api.o
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) -shared -o $@ $^ $(LDFLAGS)

clean:
	$(RM) -rf build lib

clean_all: clean

# install python-package
Pyinstall: ${PSM_DYLIB}
	rm -rf python-package/pyprimal/lib/
	mkdir python-package/pyprimal/lib/
	cp -rf ${PSM_DYLIB} python-package/pyprimal/lib/
	cd python-package; python setup.py install; cd ..

# Script to make a clean installable R package.
Rpack:
	$(MAKE) clean
	rm -rf primal PRIMAL_1.0.tar.gz
	cp -r R-package primal
	rm -rf primal/src/*.o primal/src/*.so primal/src/*.dll
	rm -rf primal/src/*/*.o
	rm -rf primal/demo/*.model primal/demo/*.buffer primal/demo/*.txt
	rm -rf primal/demo/runall.R
	cp -r src primal/src/src
	cp -r include primal/src/include
	cat R-package/src/Makevars.in|sed '2s/.*/PKGROOT=./' | sed '3s/.*/ENABLE_STD_THREAD=0/' > primal/src/Makevars.in
	cp primal/src/Makevars.in primal/src/Makevars.win
	cp primal/src/Makevars.in primal/src/Makevars

Rbuild:
	$(MAKE) Rpack
	R CMD build --no-build-vignettes primal
	rm -rf primal

Rcheck:
	$(MAKE) Rbuild
	R CMD check  PRIMAL_1.0.0.tar.gz

Rinstall:
	$(MAKE) Rbuild
	R CMD INSTALL  PRIMAL_1.0.0.tar.gz

-include build/*.d
-include build/*/*.d
