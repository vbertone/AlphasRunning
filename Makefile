default:
	pdflatex AlphasRunning.tex
	pdflatex AlphasRunning.tex
	#bibtex   AlphasRunning
	#pdflatex AlphasRunning.tex
	pdflatex notes.tex && notes.tex

clean:
	rm -rf *~ *.aux *.log *.pdf.pdf *.toc *.blg *.pdf *.out *.idx *.ind *.ilg *.dvi *.ps
