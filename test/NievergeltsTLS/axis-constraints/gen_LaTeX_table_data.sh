#!/bin/bash


# ................................................................
#          Phi     Theta  Residual    Radius    Pitch1    Sweep1
# ................................................................
#
#     125.5000  131.5000    0.0000  172.5263  686.9517   32.4521
##################################################################

for i in $(seq 3 10); do
#                                                                    atom   &   radius  &   pitch   &   sweep   &  residual  & err_rad  & err_pitch  \\
grep -A3 Phi testhelix.$i.360.out |tail -n 1|awk -v atom=$i '{printf "%i\t %s\t %.1f\t %s\t %.1f\t %s\t %.1f\t %s\t %.3f\t %s\t %.2f\t %s\t %.2f\t   %s%s \n", atom, "&", $4, "&", $5, "&", $6, "&", $3, "&", (100*(195-$4)/195), "&", (100*(600-$5)/600), "\\", "\\"  }'

done

grep -A3 Phi testhelix.10.720.out |tail -n 1|awk -v atom=$i '{printf "%i\t %s\t %.1f\t %s\t %.1f\t %s\t %.1f\t %s\t %.3f\t %s\t %.2f\t %s\t %.2f\t   %s%s \n", atom, "&", $4, "&", $5, "&", $6, "&", $3, "&", (100*(195-$4)/195), "&", (100*(600-$5)/600), "\\", "\\"  }'
grep -A3 Phi testhelix.10.1080.out |tail -n 1|awk -v atom=$i '{printf "%i\t %s\t %.1f\t %s\t %.1f\t %s\t %.1f\t %s\t %.3f\t %s\t %.2f\t %s\t %.2f\t   %s%s \n", atom, "&", $4, "&", $5, "&", $6, "&", $3, "&", (100*(195-$4)/195), "&", (100*(600-$5)/600), "\\", "\\"  }'
