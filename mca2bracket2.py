#!/usr/bin/python

import re 
import os
import sys
import math

import minilib as biolib
import mdg_dt as MDGDT

ix_pattern = re.compile('([A-Za-z]|\'[0-9]\')-?[0-9]+(\.[A-Z])?')

def is_ss_pair(line, edge_re, pair_re):
    # 'antiparallel' in line -- not in online version
    return 'cis' in line and \
        edge_re.search(line) and \
        pair_re.search(line)

###
def read_data(fn):
    sequence = []
    structure = []

    unidentified = set([])
    edge_re = re.compile('W[a-z]/W[a-z]')
    pair_re = re.compile('([AG]-U)|(C-G)|(G-[CU])|(U-[AG])')

    mode = 'ignore'
    for line in open(fn):
        line = line.strip()
        if not line: 
            mode = 'ignore'
        elif line.endswith('----'):
            if line.startswith('Residue conformations'):                
                mode = 'sequence'
                if len(sequence) > 0:
                    """ ignore multimodels """
                    break
            elif line.startswith('Base-pairs'):
                mode = 'structure'
            else:
                mode = 'ignore'
                pass
            continue
            pass        

        if mode == 'ignore': 
            continue
        elif mode == 'sequence':            
            line = line.split()
            base = line[2].strip()
            if base in biolib.BASEDICT:
                base = biolib.BASEDICT[base]
                sequence.append((line[0], base))
            else:
                unidentified.add(base)                
        elif mode == 'structure':
            # print line, mode
            # print 'edge_re:', edge_re.search(line)
            # print 'pair_re:', pair_re.search(line)

            if is_ss_pair(line, edge_re, pair_re):
                line = line.split()
                line = line[0]


                try:
                    base1 = ix_pattern.search(line)
                    line = line[base1.end():]
                    base2 = ix_pattern.search(line)
                    # print base1.group(), base2.group()
                    structure.append([base1.group(), base2.group()])
                except:
                    print line, 'fail'
                    sys.exit(1) #pass

                """
                if len(tmp) == len(line) - 1:
                    structure.append(line[0].split('-'))
                elif len(tmp) == len(line) - 2:
                    pass
                else:
                    p = line.find('-')
                """
            
        else:
            continue
        
        pass

    if len(unidentified) > 0:
        sys.stderr.write('UNIDENTIFIED RESIDUES: %s\n' % \
                         str(list(unidentified)))
        pass

    # print structure
    return sequence, structure



###
def create_brackets(structure, sequence):    
    brackets = ['.' for i in xrange(len(sequence))]

    index = [x[0] for x in sequence]
    
    structure = [(index.index(b1), index.index(b2))
                 for b1, b2 in structure]
    ss_tree = MDGDT.MDG_Stem()
    ss_tree.assemble(sorted(structure), outp=sys.stderr)
    
    knotless = sorted(ss_tree.preorder())
    if ss_tree.has_tertiary_pairs() and ss_tree.has_tertiary_stem():
        print ss_tree.te_contacts
        return None
    # print 'KLESS', knotless

    def write_structure(structure, n):
        brackets = ['.' for i in xrange(n)]
        
        for i, j in structure:
            if brackets[i] != '.' or brackets[j] != '.':
                sys.stderr.write('%s %s %i %i already occupied\n' % \
                                 (b1, b2, i, j))
                if brackets[i] == '(' or brackets[i] == '<':
                    brackets[i] = '<'
                elif brackets[i] == ')' or brackets[i] == '>':
                    brackets[i] ='#'

                if brackets[j] == '(' or brackets[j] == '<':
                    brackets[j] = '@'
                elif brackets[j] == ')' or brackets[j] == '>':
                    brackets[j] ='>'

            else:
                brackets[i] = '('
                brackets[j] = ')'
                
            pass

        # sys.stderr.write('***\n')
        return ''.join(brackets)
    
    return (write_structure(structure, len(sequence)), 
            write_structure(knotless, len(sequence)))

###
def main(argv):
#
    if len(argv) < 1:
        sys.stderr.write('Usage: python mca2bracket2.py <MC-Annotate-file>\n')
        sys.exit(1)
    elif not os.path.exists(argv[0]):
        sys.stderr.write('Error: file %s does not exist. Aborting.\n' % argv[0])
        sys.exit(1)
    
    seq, struc = read_data(argv[0])
    brackets = create_brackets(struc, seq)

    print 'Structure with knots:\n%s\nKnotless:\n%s\n' % brackets


    return None

if __name__ == '__main__': main(sys.argv[1:])
