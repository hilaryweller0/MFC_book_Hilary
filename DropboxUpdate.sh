lyx --export pdflatex -f all main.lyx; \
    pdflatex main.tex; bibtex main; pdflatex main; pdflatex main
rmtexall main.tex
cp main.pdf /home/hilary/Dropbox/MFC_book/climate_numFigures/HilaryMetNumerics.pdf

dropboxd &

#cp -u metNumerics.tex /home/hilary/Dropbox/MFC_book/climate_num.tex
meld metNumerics.tex /home/hilary/Dropbox/MFC_book/climate_num.tex &
cp -ru climate_numFigures /home/hilary/Dropbox/MFC_book

