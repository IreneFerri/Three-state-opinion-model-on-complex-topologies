# -------------------------------------------------------
clear
reset

# -------------------------------------------------------
# Titles 
# 
set title '|m| peak - {/Symbol a} = 0' font 'Times-Roman, 45'
unset title
set title  '{/Symbol a} = 0' font 'Times-Roman, 45' offset 0, -2

k=16
#unset label
set label '(c)' at 2, 0.50 font 'Times-Roman, 45'

# Margins
#set lmargin 15
#set tmargin 1

# Axis

set ylabel 'Position <|m|>_{max}' font 'Times-Roman, 45' offset -5, 0
unset ylabel
set xlabel 'k_{in}' font 'Times-Roman, 45' offset 0, -0.8
#unset xlabel
# Tics

set xtics font 'Times-Roman, 45' offset 0, -0.4 
#set xtics ('0' 0, '' 0.4, '' 0.8, '' 1.2, '' 1.6)
#set xtics add ('0.66' 0.66) scale 2
#set xtics out
set ytics 0,0.2,1 font  'Times-Roman, 45'
set y2tics 0,0.2,1 font  'Times-Roman, 45'
#set ytics 0,1,1
#unset xtics
#unset ytics

# Key

set key font 'Times-Roman, 40'
#set key inside top vertical left 
set key at 12.5, 0.9
#set key box

#unset key


# Ranges
set xrange [1:15]
set yrange [-0.05:1.05]

# Ratio and border
set border lw 2
#set size 0.80,1
set size square

#Margins
#set tmargin 2
#set rmargin 5
set lmargin 15
#set bmargin 5

set terminal aqua size 680, 700

#set label '11.33' at 11.5, 0.1 font 'Times-Roman, 18'
#set arrow from p1,0 to p1,0.6 nohead ls 0

# Plots communities  --------------------------------------------------------------------

plot 'newmann_a0_medT.csv' u ($1):($2):($4) title 'T/z(<|m|>_{max})' w e ls 6 
replot 'newmann_a0_medT.csv' u ($1):($3):(sqrt(($5)/((128*128)*16))) title '<|m|>_{max}' w e ls 7 

replot 'newmann_a0_medT.csv' u ($1):($2) notitle w l ls 6 lw 2
replot 'newmann_a0_medT.csv' u ($1):($3) notitle w l ls 7 lw 2
