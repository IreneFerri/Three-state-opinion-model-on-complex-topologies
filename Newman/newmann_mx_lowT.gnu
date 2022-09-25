# -------------------------------------------------------
clear
reset

# -------------------------------------------------------
# Titles 
# 
set title 'Newmann communities - T = 0.001 - 2000 steps - Rep.=100' font 'Times-Roman, 30'
#unset title
set title 'T = 0.001' font 'Times-Roman, 35' offset 0, -2

k=16
set label '(a)' at 2, 0.8 font 'Times-Roman, 35'

# Axis

set ylabel '<|m|>' font 'Times-Roman, 35' offset -3,0
set y2label 'Modularity' font 'Times-Roman, 35'   offset 5,0
#unset ylabel
set xlabel 'k_{in}' font 'Times-Roman, 35' offset 0, -1
#unset xlabel
# Tics

set xtics font 'Times-Roman, 35' offset 0, -0.4 
#set xtics ('0' 0, '' 0.4, '' 0.8, '' 1.2, '' 1.6)
#set xtics add ('0.66' 0.66) scale 2
#set xtics out
set ytics 0,0.2,1 font  'Times-Roman, 35'
set y2tics 0,0.2,1 font 'Times-Roman, 35'
#set ytics 0,1,1
#unset xtics
#unset ytics

# Key

set key font 'Times-Roman, 35' samplen 2
#set key below horizontal
set key vertical outside right
#set key box

#unset key

# Margin
set lmargin 13
set bmargin 5
set tmargin 3
#set rmargin 15

set terminal aqua size 1000, 650

# Ranges
#set xrange [0:2]
unset xrange
set yrange [-0.05:1.05]

# Ratio and border

#set size 1,1.1 

# Plots communities  --------------------------------------------------------------------
plot 'newmann_a0_lowT.csv' u ($14):(($9)/($1)):(sqrt(($13)/(($2)*16))) title '{/Symbol a} = 0'  w e ls 1 ps 3
#replot 'newmann_a0.50_lowT_new.csv' u ($14):(($9)/($1)):(sqrt(($13)/(($2)*16))) title '{/Symbol a} = 0.5'  w e ls 2 ps 3
#replot 'newmann_a0.5_lowT.csv' u ($14):(($9)/($1)):(sqrt(($13)/(($2)*16))) title '{/Symbol a} = 0.5'  w e ls 2
#replot 'newmann_a0.75_lowT_new.csv' u ($14):(($9)/($1)):(sqrt(($13)/(($2)*16))) title '{/Symbol a} = 0.75'  w e ls 3
#replot 'newmann_a0.85_lowT_new.csv' u ($14):(($9)/($1)):(sqrt(($13)/(($2)*16))) title '{/Symbol a} = 0.85'  w e ls 5
#replot 'newmann_a0.95_lowT_new.csv' u ($14):(($9)/($1)):(sqrt(($13)/(($2)*16))) title '{/Symbol a} = 0.95'  w e ls 8 
replot 'newmann_a1.00_lowT.csv' u ($14):(($9)/($1)):(sqrt(($13)/(($2)*16))) title '{/Symbol a} = 1'  w e ls 4 ps 3
replot 'newmann_a1.25_lowT_new.csv' u ($14):(($9)/($1)):(sqrt(($13)/(($2)*16))) title '{/Symbol a} = 1.25'  w e ls 7 ps 3
replot 'modularity_newmann.csv' title 'Modularity' w l ls 2 lw 3


