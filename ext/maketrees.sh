for f in *.out/*b; do python ext/plot_root.py ${f}.root < $f; done
