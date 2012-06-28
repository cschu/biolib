#!/usr/bin/python

import os
import sys
import math

STANDARD_BASES = ['A', 'C', 'G', 'U']

BASE_EDGE_DECODE = {0: '', 1: 'W', 2: 'H', 3: 'S',
                    4: 'B', 5: 'O', 6: 'C'}
BPAIR_CONFIG_DECODE = {0: '', 1: 'cis', 2: 'trans'}

"""
def decode_contact(contact):
    return BPAIR_CONFIG_DECODE[(contact>>1)&3] + \
           BASE_EDGE_DECODE[(contact>>6)] + \
           BASE_EDGE_DECODE[(contact>>3)&7]+ \
           str(contact & 1)
"""

###
def decode_contact(contact):
    return (BASE_EDGE_DECODE[contact >> 6],
            BASE_EDGE_DECODE[(contact >> 3) & 7],
            BPAIR_CONFIG_DECODE[(contact >> 1) & 3])


def decode_sequence(sequence):
    return ''.join([STANDARD_BASES[(x>>6)-1] for x in sequence])

def is_std_contact(contact):
    return (contact>>6) in [1,2,3] and \
           ((contact>>3)&7) in [1,2,3] and \
           ((contact>>1)&3) in [1,2]
    

###
def main(argv):
    return None

if __name__ == '__main__': main(sys.argv[1:])
