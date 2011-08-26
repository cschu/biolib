#!/usr/bin/python
#
# make_svm_data.py
#
# Functions for preparing data for libSVM experiments.
#
# Currently specific to AGO/miR experiments.
#

__author__ = 'Christian Schudoma'
__copyright__ = 'Copyright 2010-2011, Christian Schudoma'
__credits__ = []
__license__ = 'None'
__version__ = '0.1a'
__maintainer__ = 'Christian Schudoma'
__email__ = 'schudoma@mpimp-golm.mpg.de'
__status__ = 'Development'

"""
 make_svm_data (make_d) - 
 Functions for preparing data for libSVM experiments.

  Author:
     Christian Schudoma

 (c) 2010-2011

"""

import os
import sys
import math
import copy
import random

#
decode_dic = {(-1,-1,-1,-1,1): 'A',
              (-1,-1,-1,1,-1): 'C',
              (-1,-1,1,-1,-1): 'G',
              (-1,1,-1,-1,-1): 'U',
              (1,-1,-1,-1,-1): '$'}

#
encode_dic = {'A': [-1,-1,-1,-1,1],
              'C': [-1,-1,-1,1,-1],
              'G': [-1,-1,1,-1,-1],
              'U': [-1,1,-1,-1,-1],
              '$': [1,-1,-1,-1,-1]}

#
def draw_n_numbers(n, urn):
    """
    *draw_n_numbers(n, urn)*

    Randomly draws n integers from an urn. The urn is emptied in the process.

    Arguments:
       * n - the number of integers to be drawn
       * urn - the list of integers to be drawn from

    Returns:
       * list of integers
    """

    numbers = []
    while True:
        if len(urn) == 0 or len(numbers) == n: break
        p = random.randint(0, len(urn) - 1)
        numbers.append(urn[p])
        del urn[p]
    return numbers


#
def n_partition_set(n, setsize):
    """
    *n_partition_set(n, setsize)*

    Computes random n-partitions for sets of a given size.

    Arguments:
       * n - the size of the partitions
       * setsize - the size of the set

    Returns:
       * list of n-partitions (bins)
    """

    p_size = setsize / n
    urn = range(setsize)
    bins = []
    for i in xrange(n):
        bins.append(draw_n_numbers(p_size, urn))
    i = 0
    while True:
        if len(urn) == 0: break
        bins[i] += draw_n_numbers(1, urn)
        i += 1
        if i >= len(bins): i = 0
        pass
    return bins

#
def read_data(fi):
    """
    *read_data(fi)*

    Reads data from an unformatted fasta file.
    REQUIRES SEQUENCES WITH UNIQUE IDENTIFIERS!

    Arguments:
       * a file handle

    Returns:
       * a dictionary {seqid: seqdata}
    """

    data = {}
    last = None
    while True:
        line = fi.readline()
        if not line: break
        if line[0] == '>':
            last = line[1:].strip()
        else:
            data[last] = line.strip()
        pass
    
    return data

#
def prepare_data(data, dic=None, randomize=False):
    """
    *prepare_data(data, dic, randomize=False)*
    
    Preprocesses (sequence) data for libSVM experiments.
    1) Removal of sequences shared by different classes.
    2) Removal of duplicate sequences.
    3) Appending of '$' guards for shorter sequences.
    4) Binary encoding of sequences.
    5) Optional randomization of data.
    
    Arguments:
       * a data dictionary {seqid: sequence}
       * an encoding dictionary {char: binary vector}
       * a randomization flag

    Returns:
       * a dictionary {class: [sequences]}
    """

    # shared removal
    shared = {}
    for key, seq in sorted(data.items()):
        ago = float(key.split('_')[0][3])
        shared[seq] = shared.get(seq, set()).union(set([ago]))
    
    # duplicate removal, length measuring
    lengths = set([])
    ago_data = {}
    for seq, ago_set in shared.items():
        if len(ago_set) == 1: 
            ago_id = list(ago_set)[0]
            ago_data[ago_id] = ago_data.get(ago_id, []) + [seq]
            lengths.add(len(seq))
    # print ago_data

    """ ATTENTION: no seq-randomize after this point!!! """
    # encoding
    maxlen = max(lengths)
    if not dic is None:
        for key in ago_data:
            ago_data[key] = [encode(seq, dic, maxlen) for seq in ago_data[key]]
        # print ago_data

    # optional randomization
    if randomize is None:
        pass    
    elif randomize == 'class-shuffle':
        x, y = [], []
        for k, val in ago_data.items():
            y.extend([k for v in val])
            x.extend([v for v in val])
        random.shuffle(y)

        ago_data = {}
        for k, v in zip(y, x):
            ago_data[k] = ago_data.get(k, []) + [v]
        # print ago_data
    elif randomize == 'seq-shuffle':
        pass
    elif randomize == 'random-seq':
        """ Guess what I'm doing here...seq-randomizing ... in your face line 168! """
        n_seqs = sum([len(v) for v in ago_data.values()])
        seqs = n_random_sequences(n_seqs, length=len(ago_data.values()[0][0]))
        start = 0        
        for k in ago_data:
            end = start + len(ago_data[k])
            ago_data[k] = [encode(seq, dic, maxlen) for seq in seqs[start:end]]
            start = end

    return ago_data

#
def n_random_sequences(n, length=10):
    sequences = set()
    while True:
        if len(sequences) == n: break
        seq = ['ACGU'[random.randint(0, 3)] for i in xrange(length)]
        sequences.add(''.join(seq))
    return list(sequences)

#
def decode(string, dic):
    """
    *decode(string, dic)*

    Decodes a binary vector using a dictionary.

    Arguments:
       * string - a binary vector
       * a dictionary {vector: char}

    Returns:
       * the decoded string
    """

    decoded = ''
    i = 0
    step = len(dic.keys()[0])
    while i < len(string):
        # print i, i+step, string[i:i+step]
        decoded += dic[tuple(string[i:i+step])]
        i += step
    return decoded
    
#
def encode(string, dic, maxlen):
    """
    *encode(string, dic, maxlen)*

    Encodes a string as a binary vector, 
    padding shorter sequences at the 3'-end with '$'-guards.

    Arguments:
       * a string
       * a dictionary {char: vector}
       * the maxlen used to determine the number of 3'-$-guards

    Returns:
       * the encoded string
    """

    encoded = []
    for c in string:
        encoded.extend(dic[c])
    for i in xrange(maxlen - len(string)):
        encoded.extend(dic['$'])
    return encoded
    
#
def balance_data(data, minsize):     
    """
    *balance_data(data, minsize)*

    Balances a data set with subsets of unequal sizes.
    
    Arguments:
       * a data dictionary {class: [sequences]}
       * minsize - the size of the smallest class

    Returns:
       * the balanced data dictionary {class: [sequences]}
    """

    balanced = {}
    for key, val in copy.deepcopy(data.items()):
        if len(val) == minsize:
            balanced[key] = [v for v in val]
        else:
            numbers = draw_n_numbers(minsize, range(len(val)))
            balanced[key] = [val[ix] for ix in numbers]
        pass
    return balanced

#
def make_set(data, balanced_set=True, training_fraction=0.5):
    """
    *make_set(data, balanced_set=True, training_fraction=0.5)*

    Creates test and training sets from a data set.

    Arguments:
       * a data dictionary {class: [sequences]}
       * a balancing flag
       * a parameter determining the size of the training set

    Returns:
       * 4 lists: training labels/features, test labels/features
    """

    minsize = min([len(val) for key, val in data.items()])
    if balanced_set:
        dataset = balance_data(data, minsize)
    else:
        dataset = copy.deepcopy(data)
        pass   

    training_y, training_x = [], []
    test_y, test_x = [], []
    
    for k, val in dataset.items():
        training_size = int(math.ceil(len(val) * training_fraction))
        training_ = draw_n_numbers(training_size, range(len(val)))
        for i, v in enumerate(val):
            if i in training_:
                training_y.append(k)
                training_x.append(v)
            else:
                test_y.append(k)
                test_x.append(v)
        pass
    
    return training_y, training_x, test_y, test_x

#
def write_set(y, x, fo):
    """
    *write_set(y, x, fo)*
    
    Writes a data set (labels, features) to a file.

    Arguments:
       * set labels
       * set features
       * a (writing) file handle
    """
    
    for yx in zip(y, x):
        row = ['%i:%i' % (i, int(xx)) for i, xx in enumerate(yx[1])]
        fo.write('%s\n' % ' '.join([str(y)] + row))
    return None
            
  
###
def main(argv):

    data = read_data(open(argv[0]))
    data = make_set(data)

    for k, v in sorted(data.items()):
        print k, len(v)
        for seq in sorted(list(v)):
            print seq

    return None

if __name__ == '__main__': main(sys.argv[1:])


