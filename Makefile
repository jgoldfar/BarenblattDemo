LATEX=latexmk -pdf -bibtex
CHKTEX_ARGS=-n 3 -n 6
TEX_SOURCES=$(wildcard *.tex)
TEX_SOURCES_BASENAME=$(basename $(TEX_SOURCES))

%.pdf: %.tex
	$(LATEX) $<

clean-%: %.tex
	$(LATEX) -c $<

clean-all-%: %.tex
	$(LATEX) -C $<

clean-sources: $(addprefix clean-,$(TEX_SOURCES_BASENAME))

clean-all-sources: $(addprefix clean-all-,$(TEX_SOURCES_BASENAME))

clean: clean-fmt clean-check clean-sources

clean-all: clean clean-all-sources

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
