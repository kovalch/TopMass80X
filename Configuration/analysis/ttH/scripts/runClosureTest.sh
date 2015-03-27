#!/bin/sh

source $(dirname `readlink -f $0`)/parallelTools.sh

# FIXME: Need to check which steering parameters are not allowed, to avoid mis-use

for slope in 0.3 -0.3 0.6 -0.6 0.8 -0.8; do
    echo "Running samples reweighted with slope: $slope"
    echo
    for channel in ee emu mumu; do 
        for part in 101 201 102 202 103 203; do
            w
            $LA -f ttbarsignalplustau.root -p $part --reweightName 1st_add_bjet_pt --reweightSlope $slope -c $channel $@ &
            $LA -f ttbarsignalplustau.root -p $part --reweightName 1st_add_bjet_eta --reweightSlope $slope -c $channel $@ &
        done
        
        for part in 103 203; do
            w
            $LA -f ttbarsignalplustau.root -p $part --reweightName 2nd_add_bjet_pt --reweightSlope $slope -c $channel $@ &
            $LA -f ttbarsignalplustau.root -p $part --reweightName 2nd_add_bjet_eta --reweightSlope $slope -c $channel $@ &
            $LA -f ttbarsignalplustau.root -p $part --reweightName add_bjet_Mjj --reweightSlope $slope -c $channel $@ &
            $LA -f ttbarsignalplustau.root -p $part --reweightName add_bjet_dR --reweightSlope $slope -c $channel $@ &
        done
    done
    echo "Finished running samples reweighted with slope: $slope"
    echo
done

wait

echo
echo "Done running samples for closure test!"

