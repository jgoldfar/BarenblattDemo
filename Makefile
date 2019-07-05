SHELL=/bin/bash

get-tlmgr-deps:
	tlmgr init-usertree
	tlmgr update --self
	tlmgr install chktex latexmk nag natbib pgf pgfplots baskervald ec xkeyval nfssext-cfr ms xcolor svn-prov

make-conda-env:
	conda create -y --name ndce-python python=3
	source activate ndce-python && \
	pip install -r requirements.txt

check: doc/Makefile
	$(MAKE) -C $(dir $<) check l-check-compile
