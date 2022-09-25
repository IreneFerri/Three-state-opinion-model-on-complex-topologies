# -------
# Gnuplot Script for the visualization of the magnetization analyzed by 
# alpha3states.py (mean field) and and alpha3states_2D.py (lattice 2D)
# -----------------------------------------------------------------------------
reset
clear
clear
clear
clear
#
set multiplot layout 2,2 rowsfirst
set yrange [0:1]
#
#
set macros
#
# Placement of the a,b,c,d labels in the graphs
POS = "at graph 0.35,0.8 font 'Times-Ronam, 25'"
POS3 = "at graph 0.3,0.8 font 'Times-Ronam, 25'"
POS2 = "at graph 0.18,0.65 font 'Times-Ronam, 25'"

NOXTICSR = "set format x ''; unset xlabel; set xtics ('' 0, '' 0.4, '' 0.8,  '' 1.2, '' 1.6)"
NOXTICSL = "set format x ''; unset xlabel; set xtics ('' 0, '' 0.4,'' 0.8,  '' 1.2, '' 1.6)"
YTICSMX = "set ytics 0,0.2,0.8 font 'Times-Roman, 30'; set ylabel '|m|' font 'Times-Roman, 30' offset -1,0; set format y '%.1f'"
YTICSM0 = "set ytics 0,0.2,1 font 'Times-Roman, 30'; set ylabel '|m|' font 'Times-Roman, 30' offset -1,0; set format y '%.1f'"
NOYTICS = "set format y ''; unset ylabel"
XTICSL = "set xtics ('' 0, '0.4' 0.4, '' 0.8, '1.2' 1.2, '' 1.6) font 'Times-Roman, 30'; set xlabel 'T' font 'Times-Roman, 30' offset -2,0; set format x '%.1f'"
XTICSR = "set xtics  ('' 0, '0.4' 0.4, '' 0.8,  '1.2' 1.2, '' 1.6) font 'Times-Roman, 30'; set xlabel 'T' font 'Times-Roman, 30' offset 2,0; set format x '%.1f'"
#
#
# Margins for each row resp. column
#
TMARGIN = "set tmargin at screen 0.90; set bmargin at screen 0.58"
BMARGIN = "set tmargin at screen 0.58; set bmargin at screen 0.26"
LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.50"
RMARGIN = "set lmargin at screen 0.50; set rmargin at screen 0.85"
#
# Title
#
#set title 'Erdos-Renyi' font 'Times-Roman, 28' offset 0,-1
#unset key
#
# Graph a #metoo
#
set label 1 '#yotambien' @POS
set label 2 'N=2408 <k>=2.8' @POS2
@TMARGIN; @LMARGIN
@NOXTICSL; @YTICSM0
unset key
k1 = 2.8
#1
plot 'alpha_yotambien_a_0.csv' u (($3)/k1):(($9)/($1)):((sqrt($13))/(k1*($1))) title '{/Symbol a} = 0' w e ls 1 lw 2, 'alpha_yotambien_a_0.5.csv' u (($3)/k1):(($9)/($1)):((sqrt($13))/(k1*($1))) title '{/Symbol a} = 0.5' w e ls 2 lw 2, 'alpha_yotambien_a_1.csv' u (($3)/k1):(($9)/($1)):((sqrt($13))/(k1*($1))) title '{/Symbol a} = 1' w e ls 4 lw 2, 'alpha_yotambien_a_1.25.csv' u (($3)/k1):(($9)/($1)):((sqrt($13))/(k1*($1))) title '{/Symbol a} = 1.25' w e ls 7 lw 2
#
#set title 'MC phase diagram' font 'Times-Roman, 28' offset 0,-1
#
# Graph b #nochebuena
#
set xrange [0:2]
set label 1 '#nochebuena' @POS
set label 2 'N=20022 <k>=2.45' @POS2
@TMARGIN; @RMARGIN
@NOXTICSR; @NOYTICS
k2 = 2.45
#
plot 'alpha_nochebuena_a_0.dat' u (($3)/k2):(($9)/($1)):((sqrt($13))/(k2*($1))) title '{/Symbol a} = 0' w e ls 1 lw 2, 'alpha_nochebuena_a_0.5.dat' u (($3)/k2):(($9)/($1)):((sqrt($13))/(k2*($1))) title '{/Symbol a} = 0.5' w e ls 2 lw 2, 'alpha_nochebuena_a_1.dat' u (($3)/k2):(($9)/($1)):((sqrt($13))/(k2*($1))) title '{/Symbol a} = 1' w e ls 4 lw 2, 'alpha_nochebuena_a_1.25.dat' u (($3)/k2):(($9)/($1)):((sqrt($13))/(k2*($1))) title '{/Symbol a} = 1.25' w e ls 7 lw 2 
#
unset title
#
# Graph c #Marta
#
set xrange [0:2]
set label 1 '#martarovira' @POS
set label 2 'N=29110 <k>=3.67' @POS2
@BMARGIN; @LMARGIN
@XTICSL; @YTICSMX
k3 = 3.67

#set key horizontal  
#set key font "Times-Roman, 20" samplen 2 at first 1, first -0.3

#
plot 'Marta_a_0.00_new.csv' u (($3)/k3):(($9)/($1)):((sqrt($13))/(k3*($1))) title '{/Symbol a} = 0' w e ls 1 lw 2, 'Marta_a_0.50.csv' u (($3)/k3):(($9)/($1)):((sqrt($13))/(k3*($1))) title '{/Symbol a} = 0.5' w e ls 2 lw 2, 'Marta_a_1.00.csv' u (($3)/k3):(($9)/($1)):((sqrt($13))/(k3*($1))) title '{/Symbol a} = 1' w e ls 4 lw 2, 'Marta_a_1.25.csv' u (($3)/k3):(($9)/($1)):((sqrt($13))/(k3*($1))) title '{/Symbol a} = 1.25' w e ls 7 lw 2
#plot 'Marta_a_0.00.csv' u (($3)/k3):(($9)/($1)):((sqrt($13))/(k3*($1))) title '{/Symbol a} = 0' w e ls 1 lw 2,  'Marta_a_0.50.csv' u (($3)/k3):(($9)/($1)):((sqrt($13))/(k3*($1))) title '{/Symbol a} = 0.5' w e ls 2 lw 2, 'Marta_a_1.00.csv' u (($3)/k3):(($9)/($1)):((sqrt($13))/(k3*($1))) title '{/Symbol a} = 1' w e ls 4 lw 2, 'Marta_a_1.25.csv' u (($3)/k3):(($9)/($1)):((sqrt($13))/(k3*($1))) title '{/Symbol a} = 1.25' w e ls 7 lw 2
#
# Graph d #BAk2
#
set label 1 'BA tree-like' @POS
set label 2 'N=1000 <k>=2' @POS2
@BMARGIN; @RMARGIN
@XTICSR; @NOYTICS
k4 = 2
#
set xrange [0:2]
unset key
set key horizontal below 
set key box
#set key font "Times-Roman, 20" samplen 2 at second 2.9, first 1.9 
#
plot 'alpha_BA.N1000.a0_m1.dat' u (($3)/k3):(($9)/($1)):((sqrt($13))/(k4*($1))) title '{/Symbol a} = 0' w e ls 1 lw 2, 'alpha_BA.N1000.a0.5_m1.dat' u (($3)/k3):(($9)/($1)):((sqrt($13))/(k4*($1))) title '{/Symbol a} = 0.5' w e ls 2 lw 2, 'alpha_BA.N1000.a1_m1.dat' u (($3)/k3):(($9)/($1)):((sqrt($13))/(k4*($1))) title '{/Symbol a} = 1' w e ls 4 lw 2,  'alpha_BA.N1000.a1.25_m1.dat' u (($3)/k3):(($9)/($1)):((sqrt($13))/(k4*($1))) title '{/Symbol a} = 1.25' w e ls 7 lw 2
#
#
