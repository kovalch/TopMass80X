#!/bin/zsh

mkdir -p FileLists
rm FileLists/Histo*

#uncomment next line to include tt+V!! 
#foreach sample (run qcd dyee dymumu dytautau ww wz zz wtolnu ttgjets ttbarW ttbarZ single ttbarbg ttbarsignal)
foreach sample (run qcd dyee dymumu dytautau ww wz zz wtolnu single ttbarbg ttbarsignal)
   foreach channel (ee emu mumu)
     
      foreach Syst  (Nominal \
                     PERUGIA11 \
                     JES_UP JES_DOWN JER_UP JER_DOWN \ 
                     PU_UP PU_DOWN TRIG_UP TRIG_DOWN \
                     LEPT_UP LEPT_DOWN \
                     KIN_UP KIN_DOWN \
                     SCALE_UP SCALE_DOWN MATCH_UP MATCH_DOWN MASS_UP MASS_DOWN \
                     HAD_UP HAD_DOWN \
                     POWHEG POWHEGHERWIG MCATNLO \
                     BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN \
                     BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN \
                     BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN
#                      BTAG_BEFF_UP BTAG_BEFF_DOWN BTAG_CEFF_UP BTAG_CEFF_DOWN BTAG_LEFF_UP BTAG_LEFF_DOWN
)

        if [ -d selectionRoot/$Syst/$channel ] ; then
            ls -1 selectionRoot/$Syst/$channel/${channel}_$sample*.root >> FileLists/HistoFileList_$Syst\_$channel.txt
            ls -1 selectionRoot/$Syst/$channel/${channel}_$sample*.root >> FileLists/HistoFileList_$Syst\_combined.txt
        fi
      end

      foreach Sys (DY_UP DY_DOWN BG_UP BG_DOWN)
        if [ -d selectionRoot/Nominal/$channel ] ; then
            cp FileLists/HistoFileList_Nominal\_$channel.txt FileLists/HistoFileList_$Sys\_$channel.txt
            cp FileLists/HistoFileList_Nominal_combined.txt FileLists/HistoFileList_$Sys\_combined.txt
        fi
      end

   end

end
