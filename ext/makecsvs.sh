stack exec AllForTops-exe -- --outfolder data1516 --infiles data1516.infiles --lumi 1 --xsecfile data/XSection-MC15-13TeV.data > data1516/log 2>&1 &
stack exec AllForTops-exe -- --outfolder ttmc --infiles ttmc.infiles --lumi 36000 --xsecfile data/XSection-MC15-13TeV.data > ttmc/log 2>&1 &
