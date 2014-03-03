#!/bin/sh

source $(dirname `readlink -f $0`)/parallelTools.sh

for i in ee emu mumu; do 
    for s in 0.08 0.16 0.24 0.32 -0.08 -0.16 -0.24 -0.32; do
        w
        $LA -f tau.root -s closure --closure ytop --slope $s -c $i &
    done
done
w

for i in ee emu mumu; do 
    for s in 0.002 0.004 0.006 0.008 -0.002 -0.004 -0.006 -0.008; do
        w
        $LA -f tau.root -s closure --closure pttop --slope $s -c $i &
    done
done
w

for i in ee emu mumu; do
    w
    $LA -f tau.root -s closure --closure nominal -c $i &
done
w


echo "Done!"

echo "Now run unfoldClosure.sh"

