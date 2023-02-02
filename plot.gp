set autoscale

set xlabel 'Axe des x'
set ylabel 'Axe des y'

splot 'sol000.txt' w p title 'Contribution du proc 0' \
, 'sol001.txt' w p title 'Contribution du proc 1' \
, 'sol002.txt' w p title 'Contribution du proc 2' \
, 'sol003.txt' w p title 'Contribution du proc 3'

set term png
set output "numerical_solution.png"
replot
set term x11
