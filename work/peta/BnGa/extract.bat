
root << end
.L /home/ostrick/PWA/macros/Extract.cpp
.x /home/ostrick/PWA/macros/ExtractPlot.C
.q
end
pdflatex PlotMultipoles.tex

enscript -o PWAcfg.ps PWA.cfg
ps2pdf PWAcfg.ps 
pdftk PlotMultipoles.pdf PWAcfg.pdf  cat  output Result.pdf

rm PWAcfg.*

okular Result.pdf