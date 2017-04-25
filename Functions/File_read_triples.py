__author__ = 'schwahn'
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
