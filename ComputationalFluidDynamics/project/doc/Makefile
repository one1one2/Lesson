TeX=xelatex
latexfile=FVM

$(latexfile).pdf: *.tex ./images/* 
	while ($(TeX) $(latexfile); \
		grep -q "Return to get cross" $(latexfile).log ) do true;\
	done

all: $(latexfile)
$(latexfile): *.tex
	$(TeX) $(latexfile);
	bibtex $(latexfile);
	$(TeX) $(latexfile);
	$(TeX) $(latexfile);

clean:
	rm -f *.log *.out *.blg *.aux *.toc $(latexfile).dvi $(latexfile).ps
