#!/usr/bin/python

import os
import sys
import math

IMATRICES_2002 = {
    'WWc':    
    {'AA': 'I4', 'AC': 'i2', 'GU': 'i2', 'AG': 'I3',
     'CC': 'I6', 'CA': 'I2', 'CG': 'I1', 'UU': 'I6',
     'GG': None, 'GC': 'I1', 'AU': 'I1', 'GA': 'I3',
     'UG': 'I2', 'UA': 'I1', 'UC': 'I5', 'CU': 'I5'},
    'WWt':
    {'AA': 'I4', 'AC': 'I3', 'AG': None, 'CC': 'I6',
     'CA': 'I3', 'CG': 'I2', 'GU': 'I3', 'GG': 'I4',
     'UU': 'I6', 'AU': 'I1', 'UC': 'I5', 'GA': None,
     'GC': 'I2', 'UG': 'I3', 'UA': 'I1', 'CU': 'I5'},
    'WHc':
    {'AA': None, 'AC': None, 'AG': 'I3', 'CC': 'I2',
     'CA': None, 'CG': 'I1', 'GU': None, 'GG': 'I4',
     'UU': 'I2', 'AU': 'I3', 'UC': None, 'GA': 'I3',
     'GC': None, 'UG': 'I1', 'UA': 'I1', 'CU': 'I1'},
    'WHt':
    {'AA': 'I4', 'AC': None, 'AG': 'I4', 'CC': 'I1',
     'CA': 'I2', 'CG': 'I2', 'GU': 'I4', 'GG': 'I5',
     'UU': 'I2', 'AU': None, 'UC': None, 'GA': None,
     'GC': None, 'UG': 'I3', 'UA': 'I1', 'CU': None},
    'WSc':
    {'AA': 'I1', 'AC': 'I1', 'AG': 'I1', 'CC': 'I2',
     'CA': 'I2', 'CG': 'I2', 'GU': 'I3', 'GG': 'I5',
     'UU': 'I4', 'AU': 'I1', 'UC': 'I4', 'GA': 'I3',
     'GC': 'I3', 'UG': 'I4', 'UA': 'I4', 'CU': 'I2'},
    'WSt':
    {'AA': 'I1', 'AC': 'I1', 'AG': 'I1', 'CC': 'I1',
     'CA': 'I1', 'CG': 'I1', 'GU': 'I2', 'GG': None,
     'UU': 'I3', 'AU': 'I1', 'UC': 'I3', 'GA': None,
     'GC': 'I2', 'UG': 'I4', 'UA': 'I3', 'CU': 'I1'},
    'HHc':
    {'AA': None, 'AC': None, 'AG': 'I2', 'CC': None,
     'CA': None, 'CG': 'I1', 'GU': None, 'GG': 'I1',
     'UU': None, 'AU': None, 'UC': None, 'GA': 'I2',
     'GC': 'I1', 'UG': None, 'UA': None, 'CU': None},
    'HHt':
    {'AA': 'I1', 'AC': 'I1', 'AG': 'I2', 'CC': None,
     'CA': 'I1', 'CG': 'I1', 'GU': None, 'GG': 'I3',
     'UU': None, 'AU': 'I2', 'UC': 'I2', 'GA': 'I2',
     'GC': 'I1', 'UG': None, 'UA': 'I2', 'CU': 'I2'},
    'HSc':
    {'AA': 'Ix', 'AC': 'I1', 'AG': 'I1', 'CC': 'Ix',
     'CA': 'I1', 'CG': 'I1', 'GU': None, 'GG': 'I1',
     'UU': 'I1', 'AU': 'I1', 'UC': 'I1', 'GA': 'I1',
     'GC': None, 'UG': 'Ix', 'UA': 'I2', 'CU': 'Ix'},
    'HSt':
    {'AA': 'I1', 'AC': 'I1', 'AG': 'I1', 'CC': 'I1',
     'CA': 'I1', 'CG': None, 'GU': None, 'GG': 'I2',
     'UU': None, 'AU': 'I1', 'UC': None, 'GA': None,
     'GC': None, 'UG': 'I2', 'UA': 'I2', 'CU': 'I1'},
    'SSc':
    {'AA': 'I1', 'AC': 'I1', 'AG': 'I1', 'CC': 'I1',
     'CA': 'I1', 'CG': 'I1', 'GU': 'I1', 'GG': 'I1',
     'UU': 'I1', 'AU': 'I1', 'UC': 'I1', 'GA': 'I1',
     'GC': 'I1', 'UG': 'I1', 'UA': 'I1', 'CU': 'I1'},   
    'SSt':
    {'AA': 'I1', 'AC': 'I1', 'AG': 'I1', 'CC': None,
     'CA': None, 'CG': None, 'GU': 'I2', 'GG': 'I2',
     'UU': None, 'AU': 'I1', 'UC': None, 'GA': 'I2',
     'GC': 'I2', 'UG': None, 'UA': None, 'CU': None},
    'BB':
    {'AA': 'I1', 'AC': 'I1', 'AG': None, 'CC': 'I3',
     'CA': 'I2', 'CG': None, 'GU': 'I1', 'GG': 'I1',
     'UU': None, 'AU': None, 'UC': None, 'GA': None,
     'GC': None, 'UG': None, 'UA': None, 'CU': None}}

IMATRICES_2009 = {
    'WWc':
    {'AA': 'I1.4', 'AC': 'I1.2a', 'AG': 'I1.3', 'CC': 'I1.6',
     'CA': 'I1.2b', 'CG': 'I1.1', 'GU': 'I1.2a', 'GG': None,
     'UU': 'I1.7', 'AU': 'I1.1', 'UC': 'I1.5', 'GA': 'I1.3',
     'GC': 'I1.1', 'UG': 'I1.2b', 'UA': 'I1.1', 'CU': 'I1.5'},
    'WWt':
    {'AA': 'I2.7', 'AC': 'I2.4', 'AG': None, 'CC': 'I2.9',
     'CA': 'I2.3', 'CG': 'I2.6', 'GU': 'I2.4', 'GG': 'I2.7',
     'UU': 'I2.9', 'AU': 'I2.2', 'UC': 'I2.8', 'GA': None,
     'GC': 'I2.5', 'UG': 'I2.3', 'UA': 'I2.1', 'CU': 'I2.8'},
    'WHc':
    {'AA': None, 'AC': None, 'AG': 'I3.3', 'CC': 'I3.2',
     'CA': None, 'CG': 'I3.1', 'GU': None, 'GG': 'I3.4',
     'UU': 'I3.2', 'AU': 'I3.3', 'UC': None, 'GA': 'I3.3',
     'GC': None, 'UG': 'I3.1', 'UA': 'I3.1', 'CU': 'I3.2'},
    'WHt':
    {'AA': 'I4.3', 'AC': None, 'AG': 'I4.x', 'CC': 'I4.1',
     'CA': 'I4.2', 'CG': 'I4.2', 'GU': 'I4.3', 'GG': 'I4.5',
     'UU': 'I4.2', 'AU': None, 'UC': None, 'GA': None,
     'GC': None, 'UG': 'I4.4', 'UA': 'I4.1', 'CU': None},
    'WSc':
    {'AA': 'I5.1', 'AC': 'I5.1', 'AG': 'I5.1', 'CC': 'I5.2',
     'CA': 'I5.2', 'CG': 'I5.2', 'GU': 'I5.3', 'GG': 'I5.5',
     'UU': 'I5.4', 'AU': 'I5.1', 'UC': 'I5.4', 'GA': 'I5.3',
     'GC': 'I5.3', 'UG': 'I5.4', 'UA': 'I5.4', 'CU': 'I5.2'},
    'WSt':
    {'AA': 'I6.1', 'AC': 'I6.2', 'AG': 'I6.2', 'CC': 'I6.1',
     'CA': 'I6.2', 'CG': 'I6.3', 'GU': 'I6.3', 'GG': None,
     'UU': 'I6.4', 'AU': 'I6.1', 'UC': 'I6.4', 'GA': None,
     'GC': 'I6.3', 'UG': 'I6.4', 'UA': 'I6.3', 'CU': 'I6.1'},
    'HHc':
    {'AA': None, 'AC': None, 'AG': 'I7.2', 'CC': None,
     'CA': None, 'CG': 'I7.1a', 'GU': None, 'GG': 'I7.1',
     'UU': None, 'AU': None, 'UC': None, 'GA': 'I7.3',
     'GC': 'I7.1b', 'UG': None, 'UA': None, 'CU': None},
    'HHt':
    {'AA': 'I8.1', 'AC': 'I8.1', 'AG': 'I8.3', 'CC': None,
     'CA': 'I8.1', 'CG': 'I8.1', 'GU': 'I8.4', 'GG': 'I8.1',
     'UU': None, 'AU': 'I8.3', 'UC': 'I8.2', 'GA': 'I8.3',
     'GC': 'I8.2', 'UG': 'I8.2', 'UA': None, 'CU': 'I8.1'},
    'HSc':
    {'AA': 'I9.1', 'AC': 'I9.1', 'AG': 'I9.1', 'CC': 'I9.1',
     'CA': 'I9.1', 'CG': 'I9.2', 'GU': None, 'GG': 'I9.1',
     'UU': 'I9.1', 'AU': 'I9.1', 'UC': 'I9.1', 'GA': 'I9.1',
     'GC': None, 'UG': 'I9.1', 'UA': 'I9.3', 'CU': 'I9.1'},
    'HSt':
    {'AA': 'I10.1', 'AC': 'I10.1', 'AG': 'I10.1', 'CC': 'I10.1',
     'CA': 'I10.1', 'CG': None, 'GU': None, 'GG': 'I10.2',
     'UU': None, 'AU': 'I10.1', 'UC': None, 'GA': None,
     'GC': None, 'UG': 'I10.2', 'UA': 'I10.2', 'CU': 'I10.1'},
    'SSc':
    {'AA': 'I11.1', 'AC': 'I11.1', 'AG': 'I11.1', 'CC': 'I11.1',
     'CA': 'I11.1', 'CG': 'I11.1', 'GU': 'I11.1', 'GG': 'I11.1',
     'UU': 'I11.1', 'AU': 'I11.1', 'UC': 'I11.1', 'GA': 'I11.1',
     'GC': 'I11.1', 'UG': 'I11.1', 'UA': 'I11.1', 'CU': 'I11.1'},
    'SSt':
    {'AA': 'I12.1', 'AC': 'I12.1', 'AG': 'I12.1', 'CC': None,
     'CA': None, 'CG': None, 'GU': 'I12.2', 'GG': 'I12.2',
     'UU': None, 'AU': 'I12.1', 'UC': None, 'GA': 'I12.2',
     'GC': 'I12.2', 'UG': None, 'UA': None, 'CU': None}
    }


"""
Queries the isostericity matrices for all base-doublets that
are isosteric to the query doublet.
Input: iclass %c%c%c % (edge, edge, cis/trans), bases %c%c
Returns: List of isosteric doublets
"""
def get_all_isosteric_doublets(doublet, iclass, imatrix=IMATRICES_2009):

    #print doublet, iclass

    """ for A-G (WHt) special case """
    def evalIx09(ix):
        res = None
        if ix == 'I4.2' or ix == 'I4.3':
            res = ix
        return res

    def iso_equality(x, y):
        xx = x
        yy = y
        if x and x[-1] == 'x': xx = evalIx09(y)
        if y and y[-1] == 'x': yy = evalIx09(x)
        return x == y or (xx and (xx == yy))
        
    return [x for x in imatrix[iclass]
            if iso_equality(imatrix[iclass][doublet], imatrix[iclass][x])]


"""
Base pair families (i-classes) are always ordered by edge (W->H->S).
A S-H basepair has thus to be flipped.
Input: iclass %c%c%c % (edge, edge, cis/trans), bases %c%c
Returns: iclass, doublet -- possibly flipped
"""
def order_edges(doublet, iclass):

    #print doublet, iclass
    c = iclass[0], iclass[1], iclass[2][0]
    b = doublet
    flipped = False
    if (iclass[0] == 'H' and iclass[1] == 'W') or \
           (iclass[0] == 'S' and iclass[1] == 'W') or \
           (iclass[0] == 'S' and iclass[1] == 'H'):
        flipped = True
        c = iclass[1], iclass[0], iclass[2][0]
        b = doublet[1] + doublet[0]
    # print ''.join(c), b, flipped
    return ''.join(c), b, flipped

###
def main(argv):
    return None

if __name__ == '__main__': main(sys.argv[1:])
