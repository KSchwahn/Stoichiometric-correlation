__author__ = 'schwahn'
 #Python script, which reads in all the created triple files and combine it into one txt-file
 #This is necessary, if a large number of metabolites are analyzed. On most computers calculating and storing all triples in R will fill the memory and make the analysis inefficient.
import glob
import sys
filenames = glob.glob('./*T.tab')


with open(sys.argv[1], 'w') as outfile:
    fname = filenames[0]
    with open(fname, 'r') as infile:
        data = infile.read().splitlines(True)
        outfile.writelines(data)
    for fname in filenames[1:]:
        with open(fname, 'r') as infile:
            data = infile.read().splitlines(True)
            data = data[1:]
        outfile.writelines(data)
