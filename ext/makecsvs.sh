# stack exec AllForTops-exe -- --outfolder data1516 --infiles data1516.infiles --lumi 1 --xsecfile data/XSection-MC15-13TeV.data > data1516/log 2>&1 &

for f in 407342 407343 407344 410470
do
  mkdir -p $f.out
  stack exec AllForTops-exe -- --outfolder $f.out --infiles $f.infiles --lumi 36000 --xsecfile data/XSection-MC15-13TeV.data > $f.out/log 2>&1 &
done
