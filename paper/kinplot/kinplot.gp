#gnuplot file
reset

set term postscript enhanced eps color 25 dl 3
set out "kinplot.eps"

set xlabel 'x'
set ylabel 'Q^{2} (GeV2)'
set title "Kinematics of Data in NNPDFnucl1.0"
set yrange [0.9:1001]
set xrange [0.01:1.01]
set grid
set logscale x
set logscale y

set key left

plot "DATA_E139AgD.dat" u 3:4 ps 1.5 pointtype 4 t "E139AgD","DATA_NMCSnC.dat" u 3:4 ps 1.5 pointtype 4 t "NMCSnC","DATA_NMCSnCQ2.dat" u 3:4 ps 1.5 pointtype 4 t "NMCSnCQ2", "DATA_EMCSnD.dat" u 3:4 ps 1.5 pointtype 4 t "EMCSnD" 
