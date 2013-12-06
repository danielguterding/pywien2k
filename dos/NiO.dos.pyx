set terminal pdf
set output "NiO.dos.pdf"

set xticks -6,1
set xrange [-6:2]

set xlabel "energy [eV]"
set ylabel "DOS [states/eV/unit cell]"

plot "NiO.dos1ev" using 1:4 with lines title "Ni d-eg", \
     "NiO.dos1ev" using 1:5 with lines title "Ni d-t2g"