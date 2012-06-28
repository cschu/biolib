#!/usr/bin/python

import os
import sys
import math
import copy

from iso_matrix import *
from mca_coding import *


class BaseConstraint(set):
    
    def __init__(self, b1, b2, constraint):
        self.b1 = min(b1, b2)
        self.b2 = max(b1, b2)
        #self.add(constraint)
        self.constraint = constraint
        return None

    def is_clash(self, other):
        if self.b1 == other.b1:
            u, v = 0, 0
        elif self.b1 == other.b2:
            u, v = 0, 1
        elif self.b2 == other.b1:
            u, v = 1, 0
        elif self.b2 == other.b2:
            u, v = 1, 1
        else:
            return False

    pass



class SequenceModelNode(dict):

    def __init__(self):
        self.child = None
        self.parent = None
        return None
    
    pass


class IsostericSequenceGenerator():

    def __init__(self, sequence, contacts):
        self.sequence = sequence
        self.contacts = contacts
        self.results = None
        return None


    def get_sequences2(self):
        if not self.results is None:
            return self.results
        
        iso_sequences = [self.sequence]
        base_constraints = {}

        contact_matrix = {}
        contact_dic = {}
        for b1, b2, ctype in self.contacts:
            contact_matrix[b1] = contact_matrix.get(b1, []) + [b2]
            contact_matrix[b2] = contact_matrix.get(b2, []) + [b1]
            contact_dic[(b1, b2)] = ''.join(decode_contact(ctype))[:3]
            pass
        
        for b1, b2, ctype in self.contacts:
            if len(contact_matrix[b1]) > 1: print 'CLASH0', b1, b2, contact_matrix[b1]
            if len(contact_matrix[b2]) > 1: print 'CLASH1', b1, b2, contact_matrix[b2]
            pass

        sequence_model = SequenceModelNode()
        current_node = sequence_model
        model_sequence = ''
        for i, base in enumerate(self.sequence):
            current_node.child = SequenceModelNode()
            current_node = current_node.child

            if i in contact_matrix: 

                model_sequence += '-'    
                contacts = contact_matrix[i]
                if len(contacts) > 1:
                    """ In case of clash at position i """                    
                    current_node['*'] = None
                else:
                    j = contacts[0]                    
                    if len(contact_matrix[j]) > 1:
                        """ In case of clash at position j """
                        current_node['#'] = None
                    else:
                        """ No clash """
                        if i > j:
                            current_node['backtrace'] = j
                            continue
                        ctype = contact_dic[(min(i, j), max(i, j))]
                        doublet = self.sequence[i] + self.sequence[j] # -> flip?
                        iso_class, doublet, is_flipped = order_edges(doublet, ctype)
                        iso_doublets = get_all_isosteric_doublets(doublet, ctype)
                        current_node['@'] = {j: iso_doublets}
                        pass
                    pass
                

            else: 
                """ No base pairs at this position -> only one base possible (for now)"""
                model_sequence += base
                current_node[base] = None
                
                

        print model_sequence       
        current_node = sequence_model
        i = -1
        while True:
            print i, current_node.items()
            i += 1
            if current_node.child is None: break
            current_node = current_node.child
            

        return self.results


    def get_sequences(self):
        if not self.results is None:
            return self.results

        iso_sequences = [self.sequence]
        base_constraints = {}
        
        """
        For each contact, find all possible isosteric base doublets.
        For each candidate sequence from a sequence pool (starting with the source 
        sequence), generate one sequence mutant (i.e. swap the doublet with the original 
        bases) for each doublet. Then check this new candidate for clashes with 
        the sequence constraints imposed from previously established candidates.
        """
        i = 1
        for b1, b2, ctype in self.contacts:            
            # if i == 2: break
            #test_sequences = list(set([self.sequence] + iso_sequences))
            #iso_sequences = []
            test_sequences = list(set(iso_sequences))
            
            print b1, b2, i, 'of', len(self.contacts), self.sequence[b1], self.sequence[b2],
            print decode_contact(ctype)
            i += 1
            print '%i test sequences' % len(test_sequences)
            for seq in test_sequences:

                template = ''.join([seq[:b1], 
                                    '%s', 
                                    seq[b1 + 1:b2], 
                                    '%s',
                                    seq[b2 + 1:]])                                    
            
                doublet = self.sequence[b1] + self.sequence[b2]
                iso_class = ''.join(decode_contact(ctype))[:3]
                iso_class, doublet, is_flipped = order_edges(doublet, iso_class)
                iso_doublets = get_all_isosteric_doublets(doublet, iso_class)
                #print doublet, iso_class, iso_doublets
                #return None
                new_constraints = []
                for idoublet in iso_doublets:

                    if is_flipped:
                        idoublet = idoublet[::-1]
                    # print idoublet, 
                    candidate = template % (idoublet[0], idoublet[1])
                    
                    assert(len(candidate) == len(self.sequence))

                    """ 
                    Check the candidate against all constraints.
                    """
                    pattern_ok = True
                    
                    for cid, constraint in base_constraints.items():
                        if cid[0] != b1 and cid[0] != b2 and \
                                cid[1] != b1 and cid[1] != b2:
                            """ 
                            This constraint does not apply to the current pair. 
                            """
                            continue

                        """ 
                        Check current candidate sequence
                        against all base doublets of the constraint
                        until a match is found or the list is exhausted.
                        In case of a match, the candidate complies to
                        the constraint.
                        """                        
                        check_doublets = False
                        for doublet in constraint:
                            check = True
                            if cid[0] == b1:
                                #print 'X',
                                check &= candidate[b1] == doublet[0]
                            elif cid[0] == b2:
                                #print 'Y',
                                check &= candidate[b2] == doublet[1]
                            if cid[1] == b1:
                                #print 'Z',
                                check &= candidate[b1] == doublet[0]
                            elif cid[1] == b2:
                                #print 'W',
                                check &= candidate[b2] == doublet[1]
                                pass
                            if check:
                                check_doublets = True
                                break
                            pass
                        pattern_ok &= check_doublets                            
                        pass
                    
                    #print pattern_ok
                    if pattern_ok:
                        iso_sequences.append(candidate)
                        new_constraints.append((candidate[b1], candidate[b2]))
                        pass
                    pass
                key = (b1, b2, ctype)
                base_constraints[key] = base_constraints.get(key, []) + new_constraints
                pass
            pass

        self.results = sorted(list(set(iso_sequences) - set(self.sequence)))
        return self.results
                       
                                
                        
                    
"""
                    Right now, this seems like bullshit to me...
                    for cid, constraint in base_constraints.items():
                        if (cid[0] == b1 and cid[1] != b2) or \
                                (cid[1] == b2 and cid[0] != b1):
                            check_doublets = False
                            for doublet in constraint:
                                check_doublets = candidate[cid[0]] == doublet[0] and \
                                    candidate[cid[1]] == doublet[1]
                                if check_doublets: break
                                pass
                            pattern_ok &= check_doublets
                        
                        if (cid[1] == b1 and cid[0] != b2) or \
                                (cid[0] == b1 and cid[1] != b2):
                            check_doublets = False
                            for doublet in constraint:
                                check_doublets = candidate[cid[1]] == doublet[0] and \
                                    candidate[cid[0]] == doublet[1]
                                if check_doublets: break
                                pass
                            pattern_ok &= check_doublets                    
                        pass
                    """
"""                    
                    if pattern_ok:
                        contact = (b1, b2, ctype)
                        constraint = (candidate[b1], candidate[b2])

                        try:
                            base_constraints[contact].append(constraint)
                        except:
                            base_constraints[contact] = [constraint]
                            pass
                        
                        iso_sequences.append(candidate)                            
                        pass

                    pass
                pass
            pass
        
        self.results = list(set(iso_sequences) - set(self.sequence))
        return self.results

    pass
"""

###
def main(argv):
    return None

if __name__ == '__main__': main(sys.argv[1:])
