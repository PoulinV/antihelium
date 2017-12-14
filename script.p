set logscale 
set xrange [0.1:100] 
plot "build/results/pbar_TOA_spectrum.txt" u 1:2, \
"build/results/pbar_IS_spectrum.txt" u 1:2 , \
 "build/results/Debar_TOA_spectrum.txt" u 1:2, \
"build/results/Debar_IS_spectrum.txt" u 1:2 , \
 "build/results/Hebar_TOA_spectrum.txt" u 1:2, \
"build/results/Hebar_IS_spectrum.txt" u 1:2 
