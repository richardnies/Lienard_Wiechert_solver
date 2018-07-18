set term postscript
set output "plot_frontal_collision.eps"

set xlabel "x in m"
set ylabel "t in s"

set title "Frontal collision of two equal particles (q,m = 1)"

p '../output_frontal_collision.dat' every 100000 u 1:2 w lp  t 'particle 1',\
	'../output_frontal_collision.dat' every 100000 u 1:5 w lp  t 'particle 2''