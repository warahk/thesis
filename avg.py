#!/usr/bin/env python

import sys

if __name__ == '__main__':
    with open(sys.argv[1], 'r') as fp:
        print sum(float(x[x.find('[')+1:x.find(']')]) for x in fp.readlines()) / 10

    
