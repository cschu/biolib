#!/usr/bin/python
#
# mic.py
#
# Functions for mutual information computation.
#
#
#

__author__ = 'Christian Schudoma'
__copyright__ = 'Copyright 2008, Christian Schudoma'
__credits__ = []
__license__ = 'None'
__version__ = '0.1a'
__maintainer__ = 'Christian Schudoma'
__email__ = 'schudoma@mpimp-golm.mpg.de'
__status__ = 'Production - lol'

"""
 MIC - Functions for mutual information computation.

  Author:
     Christian Schudoma

 (c) 2008

"""

import os
import sys
import math


def compute_frequencies(x, pcount=0.0):
    """
    *compute_frequencies(x, pcount=0.0)*
    
    Computes the relative frequencies of elements contained in a collection.
    Optionally, a pseudocount can be specified.

    Arguments:
       * a collection
       * a non-negative pseudocount (optional; default: 0.0)

    Returns:
       * dictionary {xi: rel_freq(xi)}
    """    

    f_x = {}
    
    for xi in x:
        f_x[xi] = f_x.get(xi, int(pcount)) + 1.0
    for i in f_x:
        f_x[i] /= len(x)
    
    return f_x


def compute_joint_frequencies(x, y, pcount=0.0):
    """
    *compute_joint_frequencies(x, pcount=0.0)*
    
    Computes the relative frequencies of pairs of elements contained
    in two collections. (Calls compute_frequencies with zip(x, y).)
    Optionally, a pseudocount can be specified.

    Arguments:
       * a collection
       * a second collection
       * a non-negative pseudocount (optional; default: 0.0)

    Returns:
       * dictionary {(xi, yi): rel_freq((xi, yi))}
    """
    return compute_frequencies(zip(x, y), pcount=pcount)


def compute_mic(f_x, f_y, f_xy):
    """
    *compute_mic(f_x, f_y, f_xy)*    
    
    Computes the mutual information content between 
    two random variables from the relative frequencies of elements 
    in these variables.

    Arguments:
       * a collection of relative frequencies
       * a second collection of relative frequencies
       * a third collection of joint relative frequencies 

    Returns:
       * The mutual information using the formula:
         sum_XY[p(x,y) * log2(p(x,y)/p(x)p(y))]       
    """
    return sum([v * math.log(v / (f_x[k[0]] * f_y[k[1]]), 2.0)
                for k, v in f_xy.items()
                if f_x[k[0]] > 0.0 and f_y[k[1]] > 0.0
                and v > 0.0])

def make_alphabet(sequences):
    """
    *make_alphabet(sequences)*
    
    Computes the alphabet common to a list of sequences.

    Arguments:
       * a list of sequences

    Returns:
       * a sorted list of characters
    """
    return sorted(list(set([c for c in ''.join(sequences)])))

    



def compute_profile(sequences, alphabet=['A','C','G','U'], pcount=0.0):
    """
    *compute_profile(sequences, alphabet=['A','C','G','U'], 
                     pcount=0.0)*
    
    Computes a sequence profile for a number of sequences.
    (Calls compute_frequencies.)
    Optionally, a pseudocount can be specified.

    Arguments:
       * a list of sequences (of equal length)
       * a list of characters (an alphabet) (optional; default: RNA)
         if this is None, an alphabet is constructed from the sequences,
         otherwise, non-alphabet characters are ignored
       * a non-negative pseudocount (optional; default: 0.0)
      
    Returns:
       * list of dictionaries {xi: rel_freq(xi)}
    """
    unlimited_alphabet = False
    if alphabet is None:
        alphabet = make_alphabet(sequences)
        unlimited_alphabet = True

    maxlen = max([len(seq) for seq in sequences])
    # sequences = [seq + '-'*(maxlen-len(seq)) for seq in sequences]

    abclen = len(alphabet)
    profile = [dict(zip(alphabet, [pcount for i in xrange(abclen)]))
               for i in xrange(maxlen)]
    
    # for i in xrange(len(sequences[0])):
    for i in xrange(maxlen):
        col = [seq[i] for seq in sequences]
        if unlimited_alphabet:
            profile[i].update(compute_frequencies(col, pcount=pcount))
        else:
            col_f = compute_frequencies(col, pcount=pcount)
            known_characters = set(alphabet).intersection(set(col_f.keys()))
            for c in list(known_characters):
                profile[i][c] = col_f[c]
    
    return profile
               
    
###
def compute_mic_pos_class(sequences, classes, alphabet=['A','C','G','U']):
    """
    *compute_mic_pos_class(sequences, classes, alphabet=['A','C','G','U'])*

    Computes the mutual information content (mic) between
    each sequence column of an alignment and a class vector.

    Arguments:
       * a list of sequences (strings)
       * a list of class labels
       * an (optional) alphabet
         If this is specified, only characters in this alphabet 
         are considered for mutual information computation 
         (in order to ignore gaps).
         However, the relative column frequencies are calculated 
         based on known *and* unknown characters.

    Returns:
       * a list containing the mic for each column    
    """
    
    mic = []
    unlimited_alphabet = False
    if alphabet is None:
        alphabet = make_alphabet(sequences)
        unlimited_alphabet = True

    maxlen = max([len(seq) for seq in sequences])
    sequences = [seq + '-'*(maxlen-len(seq)) for seq in sequences]

    profile = compute_profile(sequences, 
                              alphabet=alphabet)
    f_class = compute_frequencies(classes)
    print sequences
    print classes
    print f_class
    
    for i, x in enumerate(sequences[0]):
        col_x = [seq[i] for seq in sequences]
        joint_f = compute_joint_frequencies(col_x, classes)
        if not unlimited_alphabet:
            for k in joint_f.keys():
                if not k[0] in alphabet:
                    del joint_f[k]
        print joint_f

        
        mic_i = compute_mic(profile[i], f_class, joint_f)
        mic.append(mic_i)

    return mic


"""
Example
"""

import make_svm_data as make_d
    

###
def main(argv):

    fn = argv[0]
    data = make_d.prepare_data(make_d.read_data(open(fn)))
    N = 10000

    xdata = []
    classes = []
    for x in data:
        xdata.extend([y for y in data[x]])
        classes.extend([x for y in data[x]])
    
    mic = compute_mic_pos_class(xdata, classes)
    print mic
    cmp_mic = mic
    
    rnd_mic = [0.0 for x in mic]
    for i in xrange(N):
        data = make_d.prepare_data(make_d.read_data(open(fn)),
                                   randomize='class-shuffle')
        xdata = []
        classes = []
        for x in data:
            xdata.extend([y for y in data[x]])
            classes.extend([x for y in data[x]])
        mic = compute_mic_pos_class(xdata, classes)
        rnd_mic = map(sum, zip(rnd_mic, mic))
    rnd_mic = [x/N for x in rnd_mic]

    print cmp_mic
    print rnd_mic

    fo = open('mic.csv', 'w')
    fo.write('MIC_computed,MIC_random\n')
    for xy in zip(cmp_mic, rnd_mic):
        print xy
        fo.write('%f,%f\n' % xy)
    fo.close()
    return None

if __name__ == '__main__': main(sys.argv[1:])
