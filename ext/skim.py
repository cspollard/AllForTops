from sys import argv
myargv = argv[:]
import ROOT as R

outf = R.TFile.Open(myargv[1], "CREATE")

if not outf:
  print "failed to open out file."
  exit(-1)

c = R.TChain("nominal_Loose")
c.SetBranchStatus("*_UP", 0)
c.SetBranchStatus("*_DOWN", 0)
c.SetBranchStatus("*_up", 0)
c.SetBranchStatus("*_down", 0)
c.SetBranchStatus("mc_generator_weights", 0)
cw = R.TChain("sumWeights")

for fname in myargv[2:]:
  print "adding " + fname + " to chain."
  c.Add(fname)
  cw.Add(fname)


selection = \
  "(ejets_2015_MV2c10 || ejets_2016_MV2c10 \
   || mujets_2015_MV2c10 || mujets_2016_MV2c10) \
   && (Sum$(ljet_pt > 250000) || Sum$(rcjet_pt > 250000))"


t = c.CopyTree(selection)
t.Write()

tw = cw.CopyTree("")
tw.Write()
outf.Write()

