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
    $LA -f massup -c ${chann} &
    $LA -f massdown -c ${chann} &
    $LA -f matchingup -c ${chann} &
    $LA -f matchingdown -c ${chann} &
    $LA -f scaleup -c ${chann} &
    $LA -f scaledown -c ${chann} &
    $LA -f powheg.root -c ${chann} &
    $LA -f mcatnlo -c ${chann} &
    w
done


wait

echo "Done!"

