#!/usr/bin/python
# coding: utf-8

import sys

"""
base-pair relationships
"""
""" 
default state 
"""
CT_IGN = 0 
""" 
base pairs x, y form a stack/helix 
"""
CT_STA = 1
"""
base pair x, y form a nested structure (bulge/iloop) 
"""
CT_NES = 2
""" 
base pair y follows x 
(mloop/connecting strands between stems) 
"""
CT_SEP = 3
"""
crossing base pairs (shouldn't happen 
"""
CT_CRO = 4
"""
base pair x follows y, order screwed 
(won't be dealt with)
"""
CT_PRE = 5

DEFOUT = open('/dev/null', 'w')

###
def check_relationship(contact1, contact2):
    # print contact1, contact2
    """ undecided """
    crel = CT_IGN
    sid1, sid2 = contact1[0], contact1[1]
    cid1, cid2 = contact2[0], contact2[1]
    n1, n2 = True, True
    next1, next2 = sid1 + 1, sid2 -1
        
    if next1 == cid1 and next2 == cid2:
        """ stacking/touching """
        if n1 and n2:
            crel = CT_STA
        else:
            crel = CT_NES
        pass
    elif (next1 <= cid1 and next2 > cid2) or \
            (next1 < cid1 and next2 >= cid2):
        """ separate, in same stack-sequence """            
        crel = CT_NES
    elif cid1 > sid2 and cid2 > sid2:        
        """ disjoint (separate, different stack-sequences) """        
        crel = CT_SEP
    elif cid1 < sid2 and cid2 > sid2:
        """ crossing """
        crel = CT_CRO        
    elif next1 <= sid1 - 1 and next2 <= sid2 + 1:
        """ preceeding pair """
        crel = CT_PRE 
        pass
    
    if crel == CT_IGN: print contact1, contact2
    return crel

###
class MDG_Stem :

    ###
    def __init__(self, contact=None, parent=None):

        self.contacts = []
        self.children = []
        self.ss_contacts = []
        self.te_contacts = []

        self.parent = parent        
        if contact != None :
            self.ss_contacts.append(contact)
            pass
        if parent != None :
            parent.add_child(self)
            pass
        
        return None

    ###
    def grow(self, contact, insertBack=True):
        if insertBack :
            self.ss_contacts.append(contact)            
        else :
            self.ss_contacts.insert(0, contact)
            
        return 0

    ###
    def add_child(self, stem):
        self.children.append(stem)
        return 0

    ###
    def set_parent(self, stem, override=False):
        set = False
        if self.parent == None or override:
            self.parent = stem
            set = True
            pass
        return set

    ###
    def opening_pair(self):
        pair = None
        if len(self.ss_contacts) != 0:
            pair = self.ss_contacts[0]
            pass
        return pair

    ###
    def closing_pair(self):
        pair = None
        if len(self.ss_contacts) != 0:
            pair = self.ss_contacts[-1]
            pass
        return pair

    ###
    def is_root(self):
        return self.parent == None

    ###
    def is_leaf(self) :
        return len(self.children) == 0

    ###
    def is_internal(self):
        return not self.is_leaf()

    ###
    def has_siblings(self):
        return len( self.parent.children ) > 1

    ###
    def size(self):
        return len(self.ss_contacts)

    def walk(self):
        if not self.is_root(): 
            print self.ss_contacts[0],
            print self.ss_contacts[-1]
        for child in self.children:
            child.walk()
        return None

    def preorder(self):
        pairs = []
        if not self.is_root():
            pairs.extend(self.ss_contacts)
        for child in self.children:
            pairs.extend(child.preorder())
        return pairs

    def add_tertiary_pair(self, pair):
        self.te_contacts.append(pair)
        return None

    def has_tertiary_pairs(self):
        # print 'HTE', len(self.te_contacts) > 0
        return len(self.te_contacts) > 0

    def has_tertiary_stem(self, min_stemsize=2):
        self.te_contacts.sort()
        i = 1
        stem = self.te_contacts[0:1]
        #print stem
        #print self.te_contacts
        while i < len(self.te_contacts):
            
            p1, p2 = stem[-1], self.te_contacts[i]
            #print p1, p2
            if abs(p1[0]-p2[0]) == 1 and abs(p1[1]-p2[1]) == 1:
                if len(stem) + 1 >= min_stemsize: return True
                stem.append(p2)
            else:
                stem = [p2]            
            i += 1
            pass
        return False
            
            


    ###
    def __repr__(self):
        ret = 'Dummy Stem'
        if len(self.ss_contacts) > 0:
            ret = '[ %s - %s ], %i bp\n' % \
                  (self.ss_contacts[0], self.ss_contacts[-1],
                   len(self.ss_contacts))
            pass
        return ret

    ###
    def get_stem_bases(self, contacts, max_stemsize=9999) :
        bases = []
        n = 0
        
        for c in list(contacts):
            bases.append(c[0])
            bases.append(c[1])
            n += 1
            if n == max_stemsize : break

        return bases

    ###
    def get_unpaired_bases(self, start, end, mdic) :
        bases = []
        for i in range(start, end) :
            if i in mdic : bases.append(i)
            pass
        return bases

    ###
    def assemble(self, contacts, outp=DEFOUT, min_stemsize=2):

        # print contacts
        ss_contacts = []
        nonss_contacts = []

        active_stem = self
        for c in contacts:
            # print c
            if True:
                if active_stem and active_stem.is_root():
                    # print 'XXX', c
                    stem = MDG_Stem(c, active_stem)
                    outp.write('New root-child: %s\n' % str(stem))
                    active_stem = stem
                    # print active_stem
                    ss_contacts.append(c)                    
                else:
                    acp = active_stem.closing_pair()
                    crel = check_relationship(acp, c)
                    # print 'CREL', crel, acp, c

                    if crel == CT_STA:
                        """  contact STAcks on active_stem """
                        outp.write('Growing %s\n' % active_stem)
                        active_stem.grow(c)

                    elif crel == CT_NES:
                        """ 
                        contact is a NESted pair of active_stem 
                        => new stem 
                        """                        
                        if active_stem.size() >= min_stemsize:
                            """
                            only stems >= min_stemsize are valid
                            (avoid single pairs)
                            """
                            stem = MDG_Stem(c, active_stem)
                            outp.write('Growing %s %s\n' % (stem, active_stem))
                            active_stem = stem
                            ss_contacts.append(c)
                        
                        else:
                            """ 
                            active_stem is single pair
                            traceback to parent
                            """
                            # weird format - why?
                            tmpc = [active_stem.closing_pair()]
                            stem = MDG_Stem(c, active_stem.parent)
                            # was -2 for whatever reason...
                            del active_stem.parent.children[-2]
                            active_stem = stem
                            self.te_contacts.append(tmpc[0])
                            ss_contacts.append(c)
                            pass
                        pass
                    elif crel == CT_CRO:
                         self.te_contacts.append(c)
                         pass
                    elif crel == CT_SEP:
                        """
                        contact is SEParate from active_stem
                        (e.g. multiloop siblings)
                        """
                        stem_done = False
                        if active_stem.size() < min_stemsize:
                            """
                            active_stem is single pair
                            traceback to parent
                            """
                            tmpc = [active_stem.closing_pair()]
                            par = active_stem.parent
                            # why -1 (s.above: -2)??
                            del par.children[-1]
                            self.te_contacts.append(tmpc[0])
                            active_stem = par
                            if not active_stem.is_root():
                                acp = active_stem.closing_pair()
                                crel = check_relationship(acp, c)
                               
                                if crel == CT_NES:
                                    """
                                    contact is a NESted pair of active_stem
                                    """
                                    stem_done = True
                                    stem = MDG_Stem(c, active_stem)
                                    active_stem = stem
                                    ss_contacts.append(c)
                                elif crel == CT_CRO:
                                    """
                                    CT_SEP of child can CROss parent
                                    crossing contacts => not secstr.
                                    """
                                    stem_done = True
                                    self.te_contacts.append(c)
                                    pass
                                pass
                            pass
                        outp.write('Separate\n')

                        """
                        CT_SEP of child can be separate from parent
                        => traceback
                        """
                        while not stem_done and active_stem:
                            if active_stem.is_root():
                                """ 
                                contact is direct child of root
                                stop searching
                                """
                                stem = MDG_Stem(c, active_stem)
                                outp.write('New root-child sep: %s\n' % stem)
                                active_stem = stem
                                stem_done = True
                                ss_contacts.append(c)
                            else:
                                """ continue searching """
                                active_stem = active_stem.parent
                                outp.write('Backtrace: %s\n' % active_stem)
                                if not active_stem.is_root():
                                    acp = active_stem.closing_pair()
                                    # print '----', acp, c
                                    crel = check_relationship(acp, c)
                                    if crel == CT_NES:
                                        """
                                        contact is NESted pair 
                                        of (bt) active_stem
                                        """
                                        stem = MDG_Stem(c, active_stem)
                                        outp.write('Nested (sep): %s %s\n' % \
                                                   (stem, active_stem))
                                        active_stem = stem
                                        stem_done = True
                                        ss_contacts.append(c)
                                    elif crel == CT_CRO:
                                        """ 
                                        contact CROsses (bt) active_stem
                                        """
                                        stem_done = True
                                        self.te_contacts.append(c)
                                        pass
                                    pass                                
                                else:   
                                    """
                                    we traced back to root
                                    and didn't find a fitting stem
                                    """        
                                    stem_done = True
                                    # version in mdg_dt.py
                                    # nonss_contacts.append(c)
                                    ## but shouldn't it be
                                    stem = MDG_Stem(c, active_stem)
                                    active_stem = stem
                                    ss_contacts.append(c)
                                    # ?
                                    pass
                                pass
                            pass
                        pass
                    pass
                pass            
            else:
                # --> this is never reached!!!
                """ contact is non canonical """
                self.te_contacts.append(c)
                pass
            pass # for c in contacts

        """ check if last stem is at least of min_stemsize """
        if active_stem.size() < min_stemsize:
            self.te_contacts.append(active_stem.closing_pair())
            if not active_stem.is_root():
                del active_stem.parent.children[-1]
            pass

        return 0

    
    ###
    def find_all_motifs(self, mol=None, min_stemsize=2, max_stemsize=1,
                        outp=sys.stdout):

        motifs = []
        
        if self.is_leaf() and not self.is_root():
            motifs = [['Hairpin', \
                       (self.closing_pair()[0], self.closing_pair()[1])]]
            pass
        else:
            if self.is_root():
                pass
            else:
                if len(self.children) == 1:
                    child = self.children[0]                    
                    motifs = [['Inner', (self.closing_pair()[0], \
                                         child.opening_pair()[0]), \
                               (child.opening_pair()[1], \
                                self.closing_pair()[1])]]
                    pass
                else:
                    child = self.children[0]
                    motifs = [['Multiloop',(self.closing_pair()[0], \
                                            child.opening_pair()[0])]]
                    i = 0
                    while i < len(self.children) - 1:
                        j = i + 1
                        child1 = self.children[i]
                        child2 = self.children[i+1]
                        motifs[0] += [(child1.opening_pair()[1], 
                                       child2.opening_pair()[0])]
                        i += 1
                        pass                        
                    
                    child = self.children[-1]
                    motifs[0] += [(child.opening_pair()[1], \
                                   self.closing_pair()[1])]
                    
                    pass
                pass           

            for child in self.children:
                motifs.extend(child.find_all_motifs(mol, min_stemsize, \
                                                    max_stemsize, outp))
                pass
            pass

        return motifs


def parseBracketString(s):
    stack = []
    structure = []
    for i, c in enumerate(s):
        if c == "(":
            stack.append(i)
        elif c == ")":
            try:
                structure.append((stack.pop(-1), i))
            except:
                raise ValueError("Stack prematurely empty")
    return sorted(structure)
   
    
###        
def main(argv):

    # contacts = [(0,15), (1,14), (2,13), (3,11), (4,10), (5,9)]
    # contacts = [(0,11), (1,9), (2,8), (3,7)]
    contacts = [(0,20), (1,19), (3,9), (4,8), (11,17), (12,16)]
    # 012345678901234567890
    # ((.((...)).((...)).))

    contacts = list(parseBracketString(sys.argv[1]))
    print contacts

    """ this builds the tree """
    x = MDG_Stem()
    x.assemble(contacts, outp=DEFOUT)
    for motif in x.find_all_motifs():
        print(motif)

    # print x.find_all_motifs()

    """ 
    see walk() for in-order tree-traversal 
    each node is of type MDG_Stem and contains a list ss_contacts of base pairs
    [(x,y),...]. 
    opening_pair() and closing_pair() return the first resp. last pair.
    The minimum stemsize is 2, single base pairs should be ignored.
    
    """
    # x.walk() 
    #print
    #for y in x.preorder():
    #    print y
    
    return 0

if __name__ == '__main__': main(sys.argv[1:])
