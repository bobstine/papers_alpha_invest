################################################################################
#
#
#                          GENERATING OUTPUT from Rnw 
#
#
################################################################################

RNW_SOURCES = $(wildcard *.Rnw)

LATEX_SOURCES = $(wildcard *.tex)

all: all_pdf all_R

all_html: $(RNW_SOURCES:.Rnw=.tex.html)

all_pdf: $(RNW_SOURCES:.Rnw=.pdf) 

all_pdf: $(LATEX_SOURCES:.tex=.pdf) 

all_R: $(RNW_SOURCES:.Rnw=.R)

%.tex.html: %.tex
	touch $@
	rm $@
	cat class_header $^ class_footer > $@

%.tex: %.Rnw .figures
	R CMD Sweave $<

%.R: %.Rnw
	R CMD Stangle $<

%.pdf: %.tex
	pdflatex $<

.figures:
	mkdir .figures
