# Minimal configuration for Sphinx documentation

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Dynamically export PYTHONPATH pointing to project root
PACKAGE_ROOT := $(abspath ..)
export PYTHONPATH := $(PACKAGE_ROOT):$(PYTHONPATH)

# Declare phony targets to prevent file/rule collisions
.PHONY: help Makefile clean

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Explicit clean rule to prevent Sphinx catch-all interference
clean:
	rm -rf $(BUILDDIR)/*

# Catch-all target: route all unknown targets to sphinx-build
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
