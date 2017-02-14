#!/usr/bin/env python

import os, sys
from shutil import rmtree
from subprocess import check_output

if __name__ == "__main__":
    if os.path.exists('./results'):
        rmtree('results')        
    os.mkdir('results')
    cmd = 'objs/spheres'
    methods = { 
                'serial' : './cxx_serial_spheres/%s' % cmd,
                'ispc' : './ispc_spheres/%s' % cmd,
                'cuda' : './cuda_spheres/%s' % cmd 
              }  
    resolutions = [ '256', '512', '1024' ]
    sphere_counts = [ 8, 64, 216, 512, 1000 ]

    for m in methods.keys():
        os.mkdir('results/%s' %m)
        for res in resolutions:
            for count in sphere_counts:
                sys.stdout.write("Running %s at %sx%s with %s spheres\n" % 
                                                           (m, res, res, count))
                cmd = '%s %s %s %s n' % (methods[m], res, res, count)
                outfile = 'results/%s/%sx%s_%s.txt' % (m, res, res, count)
                with open(outfile, 'w') as fp:
                    for i in range(10):
                        r = check_output(cmd, shell=True)
                        fp.write(r.strip() + '\n')

