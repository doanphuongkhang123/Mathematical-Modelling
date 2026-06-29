JOB := chemostat_food_chain_presentation
SOURCES := main.tex metadata.tex preamble.tex $(wildcard sections/*.tex)

.PHONY: all clean

all: $(JOB).pdf

$(JOB).pdf: $(SOURCES)
	latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=$(JOB) main.tex

clean:
	latexmk -C -jobname=$(JOB) main.tex
