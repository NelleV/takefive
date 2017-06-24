TEX = $(wildcard *.tex)
PDF = $(TEX:.tex=.pdf)
JUNK = $(wildcard *.log *.aux *.pdf *.blg *.bbl *.nav)


all: $(FIGURES) $(SUPPFIG) $(LATEX) $(PDF)

figures: $(FIGURES)

%.pdf: %.tex $(FIGURES) $(LATEX)
	pdflatex $*.tex
	if ( grep -q citation $*.aux ) ; then \
		bibtex $* ;\
		pdflatex $*.tex ; \
	fi
	pdflatex $*.tex

