TEX = pdflatex
TARGET = main

.PHONY: all clean

all:
	$(TEX) -interaction=nonstopmode $(TARGET).tex
	$(TEX) -interaction=nonstopmode $(TARGET).tex

clean:
	rm -f $(TARGET).aux $(TARGET).log $(TARGET).out $(TARGET).toc
