
SHELL = /bin/sh

TARGET = uq_surr.pdf
TEX    = $(TARGET:.pdf=.tex)
FIGS    = 
default: pdf
all: pdf

.IGNORE:

pdf: $(TARGET)

clean:
	/bin/rm -f *.dvi *.log *.ps *.bbl *.aux *.blg *.lof *.fff *.toc core 
	/bin/rm -f *.nav *.snm *.out *.vrb

$(TARGET): $(TEX) $(FIGS)
	pdflatex $(TEX)
	pdflatex $(TEX)
	pdflatex $(TEX)
	cp $(TARGET) ..

