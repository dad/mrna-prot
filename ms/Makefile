#BIBTEX = bibtex
#PDFLATEX = pdflatex
BIBTEX = "c:/Program Files/MiKTeX 2.9/miktex/bin/x64/bibtex.exe"
PDFLATEX = "c:/Program Files/MiKTeX 2.9/miktex/bin/x64/pdflatex.exe"
INK = "c:/Program Files (x86)/Inkscape/inkscape.exe"
INK_EPS_FLAGS = --verb EditSelectAll --verb ObjectToPath --export-eps
INK_PDF_FLAGS = --verb EditSelectAll --verb ObjectToPath --export-pdf
INK_PNG_FLAGS = --verb EditSelectAll --verb ObjectToPath --export-dpi 300 --export-background white --export-png 
#INK = "c:/Program\ Files\ \(x86\)/Inkscape/inkscape.exe"
#CONVERT = "C:/Program Files/ImageMagick-6.9.0-Q16/convert.exe"
CONVERT = convert
CONVERT_FLAGS = -flatten -density 300x300 -depth 8 -units pixelsperinch -alpha remove -compress lzw

FIG = ../figures
FINAL = ../figures

all: ms supp
bib: bibms

ms:
	$(PDFLATEX) spearman-plos.tex
	$(PDFLATEX) spearman-plos-inlinebib.tex

supp:
	$(PDFLATEX) spearman-plos-supporting.tex

rev:
#	$(PDFLATEX) response-to-reviews.tex
	$(PDFLATEX) response-to-reviews-v2.tex

bibms:
	$(PDFLATEX) spearman-plos.tex
	$(BIBTEX) spearman-plos
	$(PDFLATEX) spearman-plos.tex


finalfig: topng totiff topdf topdfsupp

topng:
	$(INK) -f $(FIG)/fig-correlations-overview.svg $(INK_PNG_FLAGS) $(FINAL)/fig1.png
	$(INK) -f $(FIG)/fig-mrna-protein-correlation.svg $(INK_PNG_FLAGS) $(FINAL)/fig2.png
	$(INK) -f $(FIG)/fig-spearman-vs-scm-sim.svg $(INK_PNG_FLAGS) $(FINAL)/fig3.png
	$(INK) -f $(FIG)/fig-nmar-modeling.svg $(INK_PNG_FLAGS) $(FINAL)/fig4.png
	$(INK) -f $(FIG)/fig-scm-vs-absolute.svg $(INK_PNG_FLAGS) $(FINAL)/fig5.png
	$(INK) -f $(FIG)/fig-exponent-final-v3.svg $(INK_PNG_FLAGS) $(FINAL)/fig6.png
	$(INK) -f $(FIG)/fig-toy-model.svg $(INK_PNG_FLAGS) $(FINAL)/fig7.png
#	$(INK) -f $(FIG)/fig_relative_dynamic_range_supp.svg $(FINAL)/figS3.png

totiff:
	$(CONVERT) $(FINAL)/fig1.png $(CONVERT_FLAGS) $(FINAL)/fig1.tiff
	$(CONVERT) $(FINAL)/fig2.png $(CONVERT_FLAGS) $(FINAL)/fig2.tiff
	$(CONVERT) $(FINAL)/fig3.png $(CONVERT_FLAGS) $(FINAL)/fig3.tiff
	$(CONVERT) $(FINAL)/fig4.png $(CONVERT_FLAGS) $(FINAL)/fig4.tiff
	$(CONVERT) $(FINAL)/fig5.png $(CONVERT_FLAGS) $(FINAL)/fig5.tiff
	$(CONVERT) $(FINAL)/fig6.png $(CONVERT_FLAGS) $(FINAL)/fig6.tiff
	$(CONVERT) $(FINAL)/fig7.png $(CONVERT_FLAGS) $(FINAL)/fig7.tiff

topdf:
	$(INK) -f $(FIG)/fig-correlations-overview.svg $(INK_PDF_FLAGS) $(FINAL)/fig-correlations-overview.pdf
	$(INK) -f $(FIG)/fig-mrna-protein-correlation.svg $(INK_PDF_FLAGS) $(FINAL)/fig-mrna-protein-correlation.pdf
	$(INK) -f $(FIG)/fig-spearman-vs-scm-sim.svg $(INK_PDF_FLAGS) $(FINAL)/fig-spearman-vs-scm-sim.pdf
	$(INK) -f $(FIG)/fig-nmar-modeling.svg $(INK_PDF_FLAGS) $(FINAL)/fig-nmar-modeling.pdf
	$(INK) -f $(FIG)/fig-scm-vs-absolute.svg $(INK_PDF_FLAGS) $(FINAL)/fig-scm-vs-absolute.pdf
	$(INK) -f $(FIG)/fig-exponent-final-v3.svg $(INK_PDF_FLAGS) $(FINAL)/fig-exponent-final.pdf
	$(INK) -f $(FIG)/fig-toy-model.svg $(INK_PDF_FLAGS) $(FINAL)/fig-toy-model.pdf

toeps:
	$(INK) -f $(FIG)/fig-correlations-overview.svg $(INK_EPS_FLAGS) $(FINAL)/fig-correlations-overview.eps
	$(INK) -f $(FIG)/fig-mrna-protein-correlation.svg $(INK_EPS_FLAGS) $(FINAL)/fig-mrna-protein-correlation.eps
	$(INK) -f $(FIG)/fig-spearman-vs-scm-sim.svg $(INK_EPS_FLAGS) $(FINAL)/fig-spearman-vs-scm-sim.eps
	$(INK) -f $(FIG)/fig-nmar-modeling.svg $(INK_EPS_FLAGS) $(FINAL)/fig-nmar-modeling.eps
	$(INK) -f $(FIG)/fig-scm-vs-absolute.svg $(INK_EPS_FLAGS) $(FINAL)/fig-scm-vs-absolute.eps
	$(INK) -f $(FIG)/fig-exponent-final-v3.svg $(INK_EPS_FLAGS) $(FINAL)/fig-exponent-final.eps
	$(INK) -f $(FIG)/fig-toy-model.svg $(INK_EPS_FLAGS) $(FINAL)/fig-toy-model.eps

topdfsupp:
	$(INK) -f $(FIG)/figS1-correlations-scatter-s-final.svg --verb EditSelectAll --verb ObjectToPath \
		--export-pdf $(FINAL)/figS1-correlations-scatter-s-final.pdf
	$(INK) -f $(FIG)/fig-exponent-supp-final.svg --verb EditSelectAll --verb ObjectToPath \
		--export-pdf $(FINAL)/figS2-exponent-supp-final.pdf
	$(INK) -f $(FIG)/fig_relative_dynamic_range_supp.svg --verb EditSelectAll --verb ObjectToPath \
		--export-pdf $(FINAL)/figS3_relative_dynamic_range_supp.pdf
	$(INK) -f $(FIG)/fig_new_vs_old_mrna.svg --verb EditSelectAll --verb ObjectToPath \
		--export-pdf $(FINAL)/figS4_new_vs_old_mrna.pdf
