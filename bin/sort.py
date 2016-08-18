#!/usr/bin/env python3.1
# Note: requires newer python

import sys

#Remove this line:
sys.argv = ('', 'file.txt')

assert(len(sys.argv) == 2)

with open(sys.argv[1], 'r') as fin:
    lines = fin.readlines()

lines_sorted = sorted(lines, key=lambda x: float(x))

back_to_file = False # Change this if needed

if back_to_file:
    with open(sys.argv[1], 'w') as fout:
        fout.writelines(lines_sorted)
else:
    for lns in lines_sorted:
        print(lns, end='') # Suppress new line
