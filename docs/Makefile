# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build
SHARDS        ?= 4
SHARD 				?= 0

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
html: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

html-nogallery:
	@echo "Not executing examples"
	@$(SPHINXBUILD) -b html -D plot_gallery=0 "$(SOURCEDIR)" "$(BUILDDIR)"/html $(SPHINXOPTS) $(0)

checksum-gallery:
	@find ../examples -type f -exec md5sum {} \; | sort -k 2  | md5sum

list-gallery:
	find ../examples -type f -name 'example*.py' | grep -v long | sort -k 2 | xargs -L 1 basename

sharded-gallery:
	$(eval TMP := $(shell mktemp -d))
	make -s list-gallery > $(TMP)/examples
	split -n l/$(SHARDS) -d $(TMP)/examples $(TMP)/examples_

	@echo "Examples in this shard:"
	$(eval SHRD := $(shell printf "$(TMP)/examples_%02d" $(SHARD)))
	@cat $(SHRD)

	@cat $(SHRD) | sed -e 's|.py|.rst|' | xargs -L 1 make

	rm source/gallery/gallery_python.zip
	rm source/gallery/gallery_jupyter.zip
	rm source/gallery/index.rst
	python3 rm-all-except.py $(SHRD)

	@rm $(TMP)/examples*
	@rmdir $(TMP)

%.rst: ../examples/%.py
	$(eval F := $(shell basename $?))
	@echo "Building example: $(F)"
	@mkdir -p source/gallery
	@touch "source/gallery/$@"
	sphinx-build $(SPHINXOPTS) -D "sphinx_gallery_conf.filename_pattern=$(F)" -D extensions=sphinx_gallery.gen_gallery -b html source build/html "source/gallery/$@"

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

