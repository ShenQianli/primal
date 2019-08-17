
clean:
	$(RM) -rf build lib bin *~ */*~ */*/*~ */*/*/*~ */*.o */*/*.o */*/*/*.o primal


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
