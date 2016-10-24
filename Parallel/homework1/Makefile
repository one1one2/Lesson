Tex=xelatex
latexfile=Gmsh

$(latexfile).pdf: *.tex
	while ($(TeX) $(latexfile); \
		grep -q "Return to get cross" $(latexfile).log ) do true;\
	done

all: $(latexfile)
$(latexfile): *.tex 
	$(TeX) $(latexfile);
	$(TeX) $(latexfile);
	$(TeX) $(latexfile);

clean:
	rm -f *.log *.out *.blg *.aux *.toc $(latexfile).dvi $(latexfile).ps
