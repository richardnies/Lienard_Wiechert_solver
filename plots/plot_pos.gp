set term postscript
set output "plot_pos.eps"

set xlabel "x in m"
set ylabel "y in m"

p '../output_pos.txt' every 10000 u 2:3 w lp lw 3 t 'particle 1',\
	'../output_pos.txt' every 10000 u 5:6 w lp lw 3 t 'particle 2''