#!/bin/zsh
#this must be run inside of Plots/Nominal
for sample in ee emu mumu combined
do
  for step in 5 6 7 8 9
  do
    #data
    tempdata=`grep Data $sample/Events$step.txt | awk -F : '{print $2;}'`
    let "data$sample$step=$tempdata"
    #tt+V boson
    tempttV=`grep "t#bar{t}+Z/W/#gamma" $sample/Events$step.txt | awk -F : '{print $2;}'`
    let "ttV$sample$step=$tempttV"
    #DY via tau
    temptau=`grep Z $sample/Events$step.txt |grep \#tau | awk -F : '{print $2;}'`
    let "tau$sample$step=$temptau"
    #DY
    tempmu=`grep Z $sample/Events$step.txt |grep \#mu | awk -F : '{print $2;}'`
    let "mu$sample$step=$tempmu"
    #diboson
    tempVV=`grep Diboson $sample/Events$step.txt | awk -F : '{print $2;}'`
    let "VV$sample$step=$tempVV"
    #wtolnu
    tempnu=`grep W+ $sample/Events$step.txt | awk -F : '{print $2;}'`
    let "nu$sample$step=$tempnu"
    #single top
    temptW=`grep Single $sample/Events$step.txt | awk -F : '{print $2;}'`
    let "tW$sample$step=$temptW"
    #ttbar bg
    tempttbg=`grep Other $sample/Events$step.txt | awk -F : '{print $2;}'`
    let "ttbg$sample$step=$tempttbg"
    # ttbar signal
    tempttsig=`grep Signal $sample/Events$step.txt | awk -F : '{print $2;}'`
    let "ttsig$sample$step=$tempttsig"
    # QCD multijet
    tempQCD=`grep QCD $sample/Events$step.txt | awk -F : '{print $2;}'`
    let "QCD$sample$step=$tempQCD"
    #Total Bg
    temptotal=`grep Total $sample/Events$step.txt | awk -F : '{print $2;}'`
    let "bgtotal$sample$step=$temptotal"
  done
done
echo "\\\begin{table}"
echo "    \\\centering \\\caption{\\\label{t-cutflow}\\\small Number of expected signal and background events, compared to the event yields in the \\\lumivalue data, after each selection step.} {\\\small"
echo "        \\\begin{tabular}{lrrrrrrr}"
echo "\\\hline $\\\mu^+\\\mu^-$ sample & 2 leptons & 2 jets & \\\met & b-tag & kin. fit \\\\\\"
echo "\\\hline"
echo "           \\\ttbar\\\ signal &     $ttsigmumu5  &       $ttsigmumu6   &       $ttsigmumu7   &       $ttsigmumu8   &       $ttsigmumu9   \\\\\\"
echo "           \\\ttbar\\\ other  &     $ttbgmumu5   &       $ttbgmumu6    &       $ttbgmumu7    &       $ttbgmumu8    &       $ttbgmumu9   \\\\\\"
echo "          $\\\ttbar+Z/W/\gamma$&     $ttVmumu5  &       $ttVmumu6   &       $ttVmumu7   &       $ttVmumu8  & $ttVmumu9  \\\\\\"
echo "           tW                 &     $tWmumu5   &       $tWmumu6    &       $tWmumu7    &       $tWmumu8    &       $tWmumu9   \\\\\\"
echo "           diboson            &     $VVmumu5   &       $VVmumu6    &       $VVmumu7    &       $VVmumu8    &       $VVmumu9   \\\\\\"
echo "           W                  &     $numumu5   &       $numumu6    &       $numumu7    &       $numumu8    &       $numumu9   \\\\\\"
echo "           DY$\\\rightarrow\\\tau\\\tau$ &     $taumumu5   &       $taumumu6    &       $taumumu7    &       $taumumu8    &       $taumumu9   \\\\\\"
echo "           DY$\\\rightarrow\\\ell^+\\\ell^-$ &     $mumumu5   &       $mumumu6    &       $mumumu7    &       $mumumu8    &       $mumumu9   \\\\\\"
echo "           QCD                  &   $QCDmumu5   &       $QCDmumu6    &       $QCDmumu7    &       $QCDmumu8    &      $QCDmumu9   \\\\\\"
echo "\\\hline"
echo "           Sum MC             &     $bgtotalmumu5   &       $bgtotalmumu6    &       $bgtotalmumu7    &       $bgtotalmumu8    &       $bgtotalmumu9   \\\\\\"
echo "           Data             &     $datamumu5   &       $datamumu6    &       $datamumu7    &       $datamumu8    &       $datamumu9   \\\\\\"

echo "\\\hline"
echo "\\\hline $\\\mu^{\\\pm}e^{\\\mp}$ sample & 2 leptons  & 2 jets & \\\met & b-tag & kin. fit\\\\\\"
echo "\\\hline"
echo "           \\\ttbar\\\ signal &     $ttsigemu5  &       $ttsigemu6   &       $ttsigemu7   &       $ttsigemu8   &       $ttsigemu9   \\\\\\"
echo "           \\\ttbar\\\ other  &     $ttbgemu5   &       $ttbgemu6    &       $ttbgemu7    &       $ttbgemu8    &       $ttbgemu9   \\\\\\"
echo "         $\\\ttbar+Z/W/\gamma$&     $ttVemu5  &       $ttVemu6   &       $ttVemu7   &       $ttVemu8  & $ttVemu9  \\\\\\"
echo "           tW                 &     $tWemu5   &       $tWemu6    &       $tWemu7    &       $tWemu8    &       $tWemu9   \\\\\\"
echo "           diboson            &     $VVemu5   &       $VVemu6    &       $VVemu7    &       $VVemu8    &       $VVemu9   \\\\\\"
echo "           W                  &     $nuemu5   &       $nuemu6    &       $nuemu7    &       $nuemu8    &       $nuemu9   \\\\\\"
echo "           DY$\\\rightarrow\\\tau\\\tau$ &     $tauemu5   &       $tauemu6    &       $tauemu7    &       $tauemu8    &       $tauemu9   \\\\\\"
echo "           DY$\\\rightarrow\\\ell^+\\\ell^-$ &     $muemu5   &       $muemu6    &       $muemu7    &       $muemu8    &       $muemu9   \\\\\\"
echo "           QCD                  &   $QCDemu5   &       $QCDemu6    &       $QCDemu7    &       $QCDemu8    &      $QCDemu9   \\\\\\"
echo "\\\hline"
echo "           Sum MC             &     $bgtotalemu5   &       $bgtotalemu6    &       $bgtotalemu7    &       $bgtotalemu8    &       $bgtotalemu9   \\\\\\"
echo "           Data             &     $dataemu5   &       $dataemu6    &       $dataemu7    &       $dataemu8    &       $dataemu9   \\\\\\"

echo "\\\hline"
echo "\\\hline \$e^{+}e^{-}$ sample & 2 leptons  & 2 jets & \\\met & b-tag & kin. fit\\\\\\"
echo "\\\hline"
echo "           \\\ttbar\\\ signal &     $ttsigee5  &       $ttsigee6   &       $ttsigee7   &       $ttsigee8   &       $ttsigee9   \\\\\\"
echo "           \\\ttbar\\\ other  &     $ttbgee5   &       $ttbgee6    &       $ttbgee7    &       $ttbgee8    &       $ttbgee9   \\\\\\"
echo "          $\\\ttbar+Z/W/\gamma$&     $ttVee5  &       $ttVee6   &       $ttVee7   &       $ttVee8  & $ttVee9  \\\\\\"
echo "           tW                 &     $tWee5   &       $tWee6    &       $tWee7    &       $tWee8    &       $tWee9   \\\\\\"
echo "           diboson            &     $VVee5   &       $VVee6    &       $VVee7    &       $VVee8    &       $VVee9   \\\\\\"
echo "           W                  &     $nuee5   &       $nuee6    &       $nuee7    &       $nuee8    &       $nuee9   \\\\\\"
echo "           DY$\\\rightarrow\\\tau\\\tau$ &     $tauee5   &       $tauee6    &       $tauee7    &       $tauee8    &       $tauee9   \\\\\\"
echo "           DY$\\\rightarrow\\\ell^+\\\ell^-$ &     $muee5   &       $muee6    &       $muee7    &       $muee8    &       $muee9   \\\\\\"
echo "           QCD                  &   $QCDee5   &       $QCDee6    &       $QCDee7    &       $QCDee8    &      $QCDee9   \\\\\\"
echo "\\\hline"
echo "           Sum MC             &     $bgtotalee5   &       $bgtotalee6    &       $bgtotalee7    &       $bgtotalee8    &       $bgtotalee9   \\\\\\"
echo "           Data             &     $dataee5   &       $dataee6    &       $dataee7    &       $dataee8    &       $dataee9   \\\\\\"

echo "\\\hline"
echo "\\\hline \$combined$ sample & 2 leptons  & 2 jets & \\\met & b-tag & kin. fit\\\\\\"
echo "\\\hline"
echo "           \\\ttbar\\\ signal &     $ttsigcombined5  &       $ttsigcombined6   &       $ttsigcombined7   &       $ttsigcombined8   &       $ttsigcombined9\\\\\\"
echo "           \\\ttbar\\\ other  &     $ttbgcombined5   &       $ttbgcombined6    &       $ttbgcombined7    &       $ttbgcombined8    &       $ttbgcombined9\\\\\\"
echo "         $\\\ttbar+Z/W/\gamma$&     $ttVcombined5  &       $ttVcombined6   &       $ttVcombined7   &       $ttVcombined8  & $ttVcombined9  \\\\\\"
echo "           tW                 &     $tWcombined5   &       $tWcombined6    &       $tWcombined7    &       $tWcombined8    &       $tWcombined9   \\\\\\"
echo "           diboson            &     $VVcombined5   &       $VVcombined6    &       $VVcombined7    &       $VVcombined8    &       $VVcombined9   \\\\\\"
echo "           W                  &     $nucombined5   &       $nucombined6    &       $nucombined7    &       $nucombined8    &       $nucombined9   \\\\\\"
echo "           DY$\\\rightarrow\\\tau\\\tau$ &     $taucombined5   &       $taucombined6    &       $taucombined7    &       $taucombined8    &       
$taucombined9\\\\\\"
echo "           DY$\\\rightarrow\\\ell^+\\\ell^-$ &     $mucombined5   &       $mucombined6    &       $mucombined7    &       $mucombined8    &       
$mucombined9\\\\\\"
echo "           QCD                  &   $QCDcombined5   &       $QCDcombined6    &       $QCDcombined7    &       $QCDcombined8    &      $QCDcombined9\\\\\\"
echo "\\\hline"
echo "           Sum MC             &     $bgtotalcombined5   &       $bgtotalcombined6    &       $bgtotalcombined7    &       $bgtotalcombined8    &       $bgtotalcombined9   \\\\\\"
echo "           Data             &     $datacombined5   &       $datacombined6    &       $datacombined7    &       $datacombined8    &       $datacombined9   \\\\\\"
echo "\\\hline"



echo "  \\\end{tabular}"
echo "  } %end small"
echo "\\\end{table}"
