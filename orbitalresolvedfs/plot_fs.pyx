set multiplot
set nodisplay
set size ratio 1.0
set width 4

set preamble r"\renewcommand{\familydefault}{\sfdefault}"

set title "GGA"

set xrange [-0.5:0.5]
set yrange [-0.5:0.5]

set xticks -0.5,0.25
set yticks -0.5,0.25

set xformat ""
set yformat ""

set grid xy

set nokey
#set key below
set keycolumns 5

set terminal pdf
set output "gga_zp.pdf"

set texthalign center
set textvalign center
set label 1 '$Z$' at 0,0
set texthalign left
set textvalign center
set label 2 '$\bar \mathsf{X}$' at 0.5,0
set label 3 '$\bar \mathsf{M}$' at 0.5,0.5

set style 17

corb1 = rgb(0.0,0.682,0.937)
corb2 = rgb(1.0,0.769,0.145)
corb3 = rgb(0.898,0.094,0.212)
corb4 = rgb(0.0,0.165,0.361)
corb5 = rgb(0.478,0.757,0.259)
set palette corb1, corb2, corb3, corb4, corb5

a=18
plot "fs.dat" select ($3==1) using 1:2 with dots color corb1 ps $4*a, \
     "fs.dat" select ($3==2) using 1:2 with dots color corb2 ps $4*a, \
     "fs.dat" select ($3==3) using 1:2 with dots color corb3 ps $4*a, \
     "fs.dat" select ($3==4) using 1:2 with dots color corb4 ps $4*a, \
     "fs.dat" select ($3==5) using 1:2 with dots color corb5 ps $4*a


set display ; refresh