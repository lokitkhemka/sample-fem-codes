# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 06:25:38 2020

@author: 91917
"""

u=8.9234555
w = 93.33
out_file = open("demo.txt","w")
out_file.write('Nodal Displacements: % 8.4f %8.4f\n' % (u,w))
out_file.write("Node    u1     u2\n")
out_file.close()