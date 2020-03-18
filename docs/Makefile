# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

html-nogallery:
	@echo "Not executing examples"
	@$(SPHINXBUILD) -b html -D plot_gallery=0 "$(SOURCEDIR)" "$(BUILDDIR)"/html $(SPHINXOPTS) $(0)

checksum-gallery:
	@find ../examples -type f -exec md5sum {} \; | sort -k 2  | md5sum

checksum-plots:
	@find source -type f -exec md5sum {} \; | sort -k 2  | md5sum

clean:
	rm -rf build/
	rm -rf source/gallery
	rm -rf source/autoapi

serve:
	@python serve.py

serve-world:
	@python serve.py --host=0.0.0.0

watch:
	@python serve.py --watch

watch-world:
	@python serve.py --host=0.0.0.0 --watch

