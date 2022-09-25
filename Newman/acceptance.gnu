# -------------------------------------------------------
clear
reset

# -------------------------------------------------------
# Titles 
# 
set title 'Newmann communities - Probability of acceptance = 0.5' font 'Times-Roman, 26'
unset title

k=16
unset label
#set label 'Marta' at 1.25, 0.75 font 'Times-Roman, 50'

# Margin
set lmargin at screen 0.15
set tmargin at screen 0.95

# Axis

set ylabel 'T' font 'Times-Roman, 50' offset -3, 0
#unset ylabel
set xlabel 'k_{in}' font 'Times-Roman, 50' offset 0, -1
#unset xlabel
# Tics

set xtics 8, 2, 16 font 'Times-Roman, 45' offset 0, -0.5 
#set xtics ('0' 0, '' 0.4, '' 0.8, '' 1.2, '' 1.6)
#set xtics add ('0.66' 0.66) scale 2
#set xtics out
set ytics font 'Times-Roman, 45'
#set ytics 0,1,1
#unset xtics
#unset ytics

# Key

set key font 'Times-Roman, 45'
set key inside bottom vertical right
#set key box

#unset key


# Ranges
set xrange [8:17]
set yrange [0:1]

# Ratio and border

set size 0.75, 1

y(x) = (16 - 2*x)/(16*log(0.5))

p1 = -(0.6*(16*log(0.5)) - 16)/2
print p1

set label '11.33' at 9.5, 0.65 font 'Times-Roman, 40'
set arrow from p1,0 to p1,0.6 nohead ls 0 lw 2

# Plots communities  --------------------------------------------------------------------

plot y(x) title 'T_{accept.}'  lw 2
replot 0.6 lw 2


