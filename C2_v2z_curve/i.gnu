set style data linespoints
set encoding iso_8859_1

set title 'PES of C_2 in cc-pVDZ basis'
#set label 'Energy extrapolation to {/Symbol e}_1={/Symbol e}_2=0' at graph .25, .95
set samples 1000
set key top right
set key spacing 1.6
set xlabel "r ({\305})"
set ylabel "Energy (Ha)"

plot \
 'energy_1_1sigma_g' u 1:2 smooth csplines lc rgb "red" title "1 ^1{/Symbol S}_g", 'energy_1_1sigma_g' u 1:2 w p lc rgb "red" notitle, \
 'energy_2_1sigma_g' u 1:2 smooth csplines lc rgb "green" title "2 ^1{/Symbol S}_g", 'energy_2_1sigma_g' u 1:2 w p lc rgb "green" notitle, \
 'energy_1_3pi_u' u 1:2 smooth csplines lc rgb "blue" title "1 ^3{/Symbol P}_u", 'energy_1_3pi_u' u 1:2 w p lc rgb "blue" notitle

set size 1.,1.; set term post eps enhanced color solid "Times-Roman" 20 ; set output 'pes_C2_2z.eps' ; replot
