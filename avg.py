#!/usr/bin/env python

import sys

with open(sys.argv[1], 'r') as fp:
    print sum(float(x.split()[5][1:-1]) for x in fp.readlines()) / 10
    
