#!/bin/bash
bed=$1
for s in seq{1..30};do
mkdir admix_${s}
cd admix_${s}
for K in 2 3 4 5; do                                                                                                                                                                                         
~/programs_that_work/admixture_linux-1.3.0/admixture --seed $s --cv ../${bed} $K | tee log${K}.out
done
grep -h CV log*.out >>../cv.out
rename admix admix_${s} *.Q
mv *.Q ../clumpak/
cd ..
done

