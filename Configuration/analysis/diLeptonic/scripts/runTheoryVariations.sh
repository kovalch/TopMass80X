#!/bin/sh

source $(dirname `readlink -f $0`)/parallelTools.sh

# #start the longest job first
# for syst in FullLeptMadgraphWithSpinCorrelation SemiLeptMadgraphWithSpinCorrelation HadronicMadgraphWithSpinCorrelation; do
#     for ch in ee emu mumu; do
#         $LA -f $syst -c $ch &
#         w
#     done
# done

# for chann in ee emu mumu; do
#     $LA -f powhegHerwig -c ${chann} -s POWHEGHERWIG &
#     w
# done

for chann in ee emu mumu; do
    $LA -f massup -s MASS_UP -c ${chann} &
    $LA -f massdown -s MASS_DOWN -c ${chann} &
    $LA -f matchingup -s MATCH_UP -c ${chann} &
    $LA -f matchingdown -s MATCH_DOWN -c ${chann} &
    $LA -f scaleup -s SCALE_UP -c ${chann} &
    $LA -f scaledown -s SCALE_DOWN -c ${chann} &
    $LA -f powheg.root -s POWHEG -c ${chann} &
    $LA -f mcatnlo -s MCATNLO -c ${chann} &
    w
done


wait

echo "Done!"

