# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS  	=
SPHINXBUILD 	= sphinx-build
SPHINXAPIDOC	= sphinx-apidoc
SPHINXPROJ    	= SuperPhot
SOURCEDIR     	= source
MODULEDIR     	= ../../PythonPackage/superphot
TESTDIR     	= ../../PythonPackage/tests
BUILDDIR    	= build
APIDOCOPTS  	= --force -M --separate
DOXYGENDIR		= ../doxygen
DOXYFILE		= Doxyfile
DOXYGEN			= doxygen

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	$(SPHINXAPIDOC) "$(MODULEDIR)" $(APIDOCOPTS) -o "$(SOURCEDIR)/_implementation" 
	$(SPHINXAPIDOC) "$(TESTDIR)" $(APIDOCOPTS) -o "$(SOURCEDIR)/_unit_tests" 
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
