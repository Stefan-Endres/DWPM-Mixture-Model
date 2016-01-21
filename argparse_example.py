#!/usr/bin/env python

from __future__ import print_function

import argparse

parser = argparse.ArgumentParser('Test argument parsing')

parser.add_argument('-c', '--compounds', nargs='+', required=True,
                    help='compounds to simulate for example'
                         'acetone water phenol')
parser.add_argument('-p', '--phases', nargs='+', required=True,
                    help='List of valid phases in equilibrium'
                         'ex. for VLE use x y')
parser.add_argument('-e', '--eos', default='DWPM',
                    choices=['DWPM'],
                    help='Equation of State / Mixture rule')
parser.add_argument('-r', type=float,
                    help='Force value of r')
parser.add_argument('-s', type=float,
                    help='Force value of s')
parser.add_argument('-m', '--model', nargs=1,
                    default="Adachi-Lu",
                    choices=['Adachi-Lu', 'Soave'],
                    help='Actvity coefficient model')
parser.add_argument('-T', '--temperature', nargs=1, type=float,
                    help='Temperature for point simulation')
parser.add_argument('-P', '--pressure', nargs=1, type=float,
                    help='Pressure for point simulation')
parser.add_argument('-z', type=float, nargs="+",
                    help='Composition for point simulation')


args = parser.parse_args()

print(args)

# Try running this with the following arguments to see the options printed out
# it's basically already in the form you are using in data_handling, but we 
# don't have to write any code.
# python argparse_example.py -c ethanol water -p x y -e DWPM

if 'T' in args:
    print("Specified point simulation")
