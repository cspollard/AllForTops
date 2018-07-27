for f in ttmc/*b; do python ext/plot_root.py $(f).root < $f; done
for f in data1516/*b; do python ext/plot_root.py $(f).root < $f; done
