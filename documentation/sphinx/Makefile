# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS  	=
SPHINXBUILD 	= sphinx-build
SPHINXAPIDOC	= sphinx-apidoc
SPHINXPROJ    	= AstroWISP
SOURCEDIR     	= source
MODULEDIR     	= ../../PythonPackage/astrowisp
TESTDIR     	= ../../PythonPackage/astrowisp/tests
BUILDDIR    	= build
APIDOCOPTS  	= --force -M --separate
DOXYGENDIR		= ../doxygen
DOXYFILE		= Doxyfile
DOXYGEN			= doxygen
INSTALLDIR      = ../../docs

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

install: html
	rm -rf $(INSTALLDIR)
	mv $(BUILDDIR)/html $(INSTALLDIR)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	./make_doc_src.sh
	mkdir -p $(DOXYGENDIR)/build/html
	cd $(DOXYGENDIR) && rm -rf build/xml && doxygen Doxyfile && cp images/* build/xml/
	rm -rf $(BUILDDIR)
	mkdir -p $(BUILDDIR)/html
	$(SPHINXAPIDOC) "$(MODULEDIR)" $(APIDOCOPTS) -o "$(SOURCEDIR)/_implementation" 
	$(SPHINXAPIDOC) "$(TESTDIR)" $(APIDOCOPTS) -o "$(SOURCEDIR)/_unit_tests" 
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	mv $(DOXYGENDIR)/build/html $(BUILDDIR)/html/doxygen
