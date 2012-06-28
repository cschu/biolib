#!/usr/bin/python

import os
import re
import sys
import glob

import minilib as biolib
import mdg_dt as MDGDT

try:
    sys.stderr = open('/tmp/sstr_error.log', 'w')
except:
    pass

ix_pattern = re.compile('([A-Za-z]|\'[0-9]\')-?[0-9]+(\.[A-Z])?')

###
def read_data(fn):
    sequence = []
    structure = []

    unidentified = set([])

    edge_re = re.compile('W[a-z]/W[a-z]')
    pair_re = re.compile('([AG]-U)|(C-G)|(G-[CU])|(U-[AG])')

    fi = open(fn)
    mode = 'ignore'
    while True:
        line = fi.readline()
        if not line: break

        line = line.strip()
        if line.endswith('----'):
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

        
        if mode == 'ignore': continue
        elif mode == 'sequence':            
            line = line.split()


            base = line[2].strip()
            if base in biolib.BASEDICT:
                base = biolib.BASEDICT[base]
                sequence.append((line[0], base))
            else:
                unidentified.add(base)
                
        elif mode == 'structure':
            is_ss_pair = 'antiparallel' in line and \
                         'cis' in line and \
                         edge_re.search(line) and \
                         pair_re.search(line)
            if is_ss_pair:
                line = line.split()
                line = line[0]


                try:
                    base1 = ix_pattern.search(line)
                    line = line[base1.end():]
                    base2 = ix_pattern.search(line)
                    structure.append([base1.group(), base2.group()])
                except:
                    print line, 'fail'
                    pass

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

    
    return sequence, structure


###
def do_something():
    return 0

###
def remove_pknots(structure):

    knotless = []
    queue = [x for x in sorted(structure)]   

    i = 0
    while True:
        if len(queue) == 0: break
        u, queue2 = queue[0], queue[1:]        
        queue = []
        for j, v in enumerate(queue2):
            if u[0] < v[1] and v[0] < u[1] and v[1] > u[1]:
                #sys.stderr.write('Pair %s is crossing pair %s\n' % \
                #                 (str(v), str(u)))
                pass
            elif u[0] == v[0] or u[1] == v[0] or \
                 u[0] == v[1] or u[1] == v[1]:
                pass
            else:
                queue.append(v)
                pass
            pass
        knotless.append(u)
        pass
    
    return knotless

###
def check_multipairs(structure):

    contacts = {}
    for p in structure:
        if p[0] in contacts:
            contacts[p[0]].append(p[1])
        else:
            contacts[p[0]] = [p[1]]
            pass
        if p[1] in contacts:
            contacts[p[1]].append(p[0])
        else:
            contacts[p[1]] = [p[0]]            
        pass

    mpairs = [(x, [str(y) for y in contacts[x]])
              for x in contacts if len(contacts[x])>1]
    for mp in mpairs:
        sys.stderr.write('%s %s\n' % (str(mp[0]), ', '.join(mp[1])))
        pass
    
    return None
    

###
def create_brackets(structure, sequence):    
    brackets = ['.' for i in xrange(len(sequence))]

    index = [x[0] for x in sequence]
    
    structure = [(index.index(b1), index.index(b2))
                 for b1, b2 in structure]
    ss_tree = MDGDT.MDG_Stem()
    ss_tree.assemble(sorted(structure), outp=open('/dev/null', 'w'))
    
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

        sys.stderr.write('***\n')
        return ''.join(brackets)

    return write_structure(structure, len(sequence)) + '\n' + \
           write_structure(knotless, len(sequence))
    

###
def main(argv):

    fo = open('darts_set_secondary_structures.fas', 'w')

    valid_ids = []
    dataset = open('../dataset')
    while True:
        line = dataset.readline()

        #print line
        if not line: break
        line = line.strip()
        if re.match('[0-9][0-9a-z]{3}(:[A-Za-z0-9]+)?', line):
            chains = ['pdb' + line[:4] + '1' + x
                      for x in line[5:]]
            if len(chains) == 1 and len(chains[0]) == 8:
                chains[0] += '?'
                pass
            # print chains
            for chain in chains:
                #print chain
                
                valid_ids.extend(glob.glob(chain + '.mca'))
            pass

        else:
            #print 'No match'
            pass

        pass

    # print valid_ids
    for pid in valid_ids:

        if not os.path.exists(pid):
            sys.stderr.write(pid + ' does not exist\n')
            continue
        sys.stderr.write(pid + '\n')
        
        seq, struct = read_data(pid)
        brackets = create_brackets(struct, seq)

        fo.write('>' + pid[:-4] + '\n')
        fo.write(''.join([x[1] for x in seq]) + '\n')
        fo.write(brackets + '\n')
        pass
    fo.close()
    
    return 0

if __name__ == '__main__': main(sys.argv[1:])




"""
1111.222....222..344455555...666....666..55555.444..31111
((((.(((....)))..(((((((((...(((....)))..))))).)))..)))))

"""
