#!/usr/bin/env python

from distutils.version import LooseVersion
import os
import subprocess
import sys

from avg import avg

if __name__ == "__main__":
    result_dir = sys.argv[1]
    output_file = sys.argv[2]
    file_names = os.listdir(result_dir + '/cuda')
    file_names.sort(key=LooseVersion)
    
    solutions = ['serial', 'ispc', 'ispc_multicore', 'cuda']  
    fp = open(output_file, 'w')
    fp.write(','.join(solutions) + '\n')
    for f in file_names:
        line = []
        line.append(f)
        for x in solutions:
            line.append(str(avg(os.path.join(result_dir, x, f))))
        fp.write(','.join(line) + '\n')
            
            
