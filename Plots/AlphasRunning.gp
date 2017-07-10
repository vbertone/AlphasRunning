 set term postscript enhanced color 'Helvetica,20' size 10,7
 set out "AlphasRunning.eps"

 set xlabel "Q^2 bin [GeV]"
 set ylabel "{/Symbol a}_s(MZ)"

 set logscale x
 set key left top

 f(x) = 0.1181
 fd(x) = 0.1170
 fu(x) = 0.1192

 set title "Simultaneous fit of {/Symbol a}_s(M_Z) and PDFs, bin by bin in Q^2"
 plot [2:1000] \
 "AlphasRunning.dat" u 1:2:3 index 0 with yerrorbars pointtype 7 pointsize 1.0 lc rgb "black" t "Fit results", \
 f(x)  w l lw 3 lc rgb "red" lt 0 t "PDG", \
 fd(x) w l lw 3 lc rgb "red" notitle, \
 fu(x) w l lw 3 lc rgb "red" notitle