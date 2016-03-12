set terminal png transparent enhanced font arial 14 size 800,600 
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb"white" behind
set grid
set logscale y
set key bottom right
set style line 3 lt 2 lc rgb "black" lw 3
set style line 2 lt 2 lc rgb "green" lw 3
set style line 1 lt 2 lc rgb "red" lw 1 
set style line 4 lt 2 lc rgb "blue" lw 3 
set style line 5 lt 3 lc rgb "orange" lw 1

set output outputname

set yrange [* : *]
set title "Shock position"
set xlabel "Time [s]"
set ylabel "Radius [cm]" 
set samples 500

p 'shockpos.txt' u 1:2 w lines ls 3 title "max v", 'shockpos.txt' u 1:4 w l ls 4 title "max div(v)"
