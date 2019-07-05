LATEX=latexmk -pdf -bibtex
PYTHON=python
PANDOC=pandoc

CHKTEX_ARGS=-n 3 -n 6
# Pick up environment's CHKTEXRC if possible
ifneq ($(CHKTEXRC),)
CHKTEX_ARGS+=--localrc=${CHKTEXRC}
endif

TEX_SOURCES=$(wildcard *.tex)
TEX_SOURCES_BASENAME=$(basename $(TEX_SOURCES))

barenblatt-profile-deps: ../requirements.txt
	pip install --user -r $<

# Generate data for images included into main file
BARENBLATT_PROFILE_SIGMA=2
images/barenblatt-profile.dat: ../barenblatt_profile.py barenblatt-profile-deps
	mkdir -p $(dir $@)
	$(PYTHON) $< $@ $(BARENBLATT_PROFILE_SIGMA) --force

images/barenblatt-profile-indep.pdf: images/barenblatt-profile-indep.tex images/barenblatt-profile.tex images/barenblatt-profile.dat

images/%-indep.tex: images/%.tex
	echo "\\documentclass{standalone}" > $@
	echo "\\usepackage{pgfplots}" >> $@
	echo "\\def\\datafiledir{.}" >> $@
	echo "\\begin{document}" >> $@
	echo "\\input{$*}" >> $@
	echo "\\end{document}" >> $@

clean-image-barenblatt-profile:
	$(RM) images/barenblatt-profile-indep.pdf images/barenblatt-profile-indep.tex
	$(RM) images/barenblatt-profile.dat

IMAGES+=barenblatt-profile
IMAGESFORMAIN+=barenblatt-profile.dat

# General target for building all image data
images: $(addprefix images/,$(IMAGESFORMAIN))

clean-images: $(addprefix clean-image-,$(IMAGES))

# Target for PDF generation (including images)
%.pdf: %.tex images
	$(LATEX) $<

%.html: %.tex images
	$(PANDOC)	$< -s --mathjax -o $@

clean-%: %.tex
	$(LATEX) -c $<

clean-all-%: %.tex
	$(LATEX) -C $<

clean-sources: $(addprefix clean-,$(TEX_SOURCES_BASENAME))

clean-all-sources: $(addprefix clean-all-,$(TEX_SOURCES_BASENAME))

clean: clean-fmt clean-check clean-sources

clean-all: clean clean-all-sources clean-images

fmt: $(addsuffix .bak,$(TEX_SOURCES_BASENAME))

clean-fmt:
	$(RM) *.bak
	$(RM) indent.log

check: $(addprefix check-,$(TEX_SOURCES_BASENAME))

l-check-compile: $(addsuffix .pdf,$(TEX_SOURCES_BASENAME))

clean-check:
	$(RM) lint-*.out

fmt-%: %.bak

%.bak: %.tex
	echo "Indenting $<"
	latexindent -w -l $<

lint-%.out: %.tex
	chktex -q $(CHKTEX_ARGS) $< 2>/dev/null | tee $@

check-%: lint-%.out
	test ! -s $<
