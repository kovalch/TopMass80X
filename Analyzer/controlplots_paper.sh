./ControlPlots -c lepton_Baseline --analysisConfig.plotsToDraw "BackgroundPlots|JESUncVariationPlots|JERUncVariationPlots|MCGeneratorPlots|PUPlots|LeptonJets|JetPts" &
./ControlPlots -c lepton_BaselinePerm --analysisConfig.normalizeToData=true -o _norm --analysisConfig.plotsToDraw="BackgroundPlots|JESUncVariationPlots|JERUncVariationPlots|MCGeneratorPlots|PUPlots|LeptonJets|BasicMasses" --analysisConfig.extraLabel "Before kinematic fit" &
./ControlPlots -c lepton --analysisConfig.normalizeToData=true -o _norm --analysisConfig.plotsToDraw="BackgroundPlots|JESUncVariationPlots|JERUncVariationPlots|MCGeneratorPlots|PUPlots|LeptonJets|BasicMasses" --analysisConfig.extraLabel "After P_{gof} selection" &
./ControlPlots -c alljets_paperHS --analysisConfig.extraLabel "After P_{gof} selection" &
