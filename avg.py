#!/usr/bin/env python

import sys

def avg(file_path):
    with open(file_path) as fp:
        x = sum(float(x[x.find('[')+1:x.find(']')]) for x in fp.readlines()) / 10
    return x

if __name__ == '__main__':
    print avg(sys.argv[1])
    
