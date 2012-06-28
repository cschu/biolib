#!/usr/bin/python

import os
import re
import sys
import glob
import pickle

# import mdg_dt2 as MDGDT

"""
Current MCA Bitformats:

pucker/seq:

XXXYYYYZZ, X: base, Y: pucker, Z: sugar configuration

stacking interactions:

XXY, X: stacking direction, Y: adjacent/non-adjacent?

contacts:

XXXYYYZZW, X: edge1, Y: edge2, Z: cis/trans conf, W: one-hbond?
"""



###
def parse_mcaindex(mca_index) :

        cid = mca_index[0]
                
        if cid == '\'':
            cid = mca_index[1]
            ix = mca_index[3:].strip()
        elif cid.isdigit() or cid == '-':
            ix = mca_index[0:].strip()
            cid = ' '
        else:
            ix = mca_index[1:].strip()
            
        ic = ' '
        if ix.find('.') != -1 :
            ic = ix[-1]
            ix = ix[:-2]
            pass
        
        return (int(ix), ic, cid)

###
def parse_indices(indices):
    pat = '([A-Za-z]|(\'[0-9]\'))?((-[0-9]{1,3})|[0-9]{1,4})(\.[A-Za-z])?'
    p = re.compile(pat)
    match = p.search(indices)
    ix1 = parse_mcaindex(match.group())
    match = p.search(indices[match.end() + 1:])
    ix2 = parse_mcaindex(match.group())

    return ix1, ix2

###
def parse_mca(line) :

    contact = {'comment': '', 'config': '', 'e1': '', 'e2': ''}

    indices, data = line.split(':')        
    contact['id1'], contact['id2'] = parse_indices(indices)
    
    data = data.split()
    for d in data:
        if d.find('-') != -1:                
            contact['b1'], contact['b2'] = d.split('-')
        elif d.find('/') != -1:
            contact['ei1'], contact['ei2'] = d.split('/')
            contact['e1'] = contact['ei1'][0]
            contact['e2'] = contact['ei2'][0]
        elif d == 'parallel' or d == 'antiparallel':
            contact['ori'] = d
        elif d == 'cis' or d == 'trans':
            contact['config'] = d
        elif (d.isupper() or d.isdigit()) and d == data[-1]:
            contact['cclass'] = d
        elif d != 'pairing':
            contact['comment'] += 'COMMENT %s\n' % d
            pass
        pass

    return contact
                
    

class MCA_Data:

    ###
    def read(self, fn):

        state_dic = {'Residue conformations': 1, \
                     'Adjacent stackings': 3, \
                     'Non-Adjacent stackings': 2, \
                     'Base-pairs': 4}

        puckers = ['unknown', \
                   'C1p_endo', 'C1p_exo', \
                   'C2p_endo', 'C2p_exo', \
                   'C3p_endo', 'C3p_exo', \
                   'C4p_endo', 'C4p_exo', \
                   'O4p_endo', 'O4p_exo', \
                   ]
        
        pucker_rdic = dict(zip(range(len(puckers)), puckers))
        pucker_dic = dict(zip(puckers, range(len(puckers))))
        

        pconf_dic = {'unknown': 0, 'anti': 1, 'syn': 2}
        pconf_rdic = {0: 'unknown', 1: 'anti', 2: 'syn'}

        base_dic = {'unknown': 0, 'A': 1, 'I': 1, 'C': 2, 'G': 3, 'T': 4, 'U': 4}
        base_rdic = {0: 'unknown', 1: 'A', 2: 'C', 3: 'G', 4: 'U'}

        stack_dic = {'upward': 0, 'downward': 1, \
                     'inward': 2, 'outward': 3}
        stack_rdic = {0: 'upward', 1: 'downward', \
                      2: 'inward', 3: 'outward'}

        edge_dic = {'': 0, 'W': 1, 'H': 2, 'S': 3, 'B': 4, 'O': 5, 'C': 6}
        edge_rdic = {0: '', 1: 'W', 2: 'H', 3: 'S', 4: 'B', 5: 'O', 6: 'C'}

        conf_dic = {'': 0, 'cis': 1, 'trans': 2}
        conf_rdic = {0: '', 1: 'cis', 2: 'trans'}

        fi = open(fn)

        state = 0
        while True:
            line = fi.readline()
            if not line: break
            # print line
            line = line.strip()

            if line.endswith('---'):
                state = state_dic[line[:line.find('--')].strip()]
            elif line.startswith('Number'):
                continue
            # pucker section            
            elif state == 1:
                index, data = line.split(':')
                index = parse_mcaindex(index.strip())
                self.resdic[index] = len(self.resdic)
                data = data.strip().split()
                # print data
                # print pucker_dic
                if not data[0] in base_dic:
                    data[0] = 'unknown'
                enc = base_dic[data[0]] << 6
                if len(data) > 2:                    
                    enc |= pucker_dic[data[1]] << 2 | pconf_dic[data[2]]
                self.sequence.append(enc)
            # stacking section
            elif state == 2 or state == 3:
                indices, data = line.split(':')
                index = tuple([self.resdic[ix] for ix in parse_indices(indices)])
                data = data.split()
                # print data
                stype = stack_dic[data[state-2].strip()] << 1 | state-2
                self.stackings[index] = stype
            # contact section
            elif state == 4:
                contact = parse_mca(line)
                #print contact,
                enc = edge_dic[contact['e1']] << 6 | \
                      edge_dic[contact['e2']] << 3 | \
                      conf_dic[contact['config']] << 1
		#print enc
                if contact['comment'].find('one_hbond') != -1:
                    enc |= 1
                    pass
                key = (self.resdic[contact['id1']], self.resdic[contact['id2']])
		print key, enc
                self.contacts[key] = enc
                pass
            pass


        fi.close()
        return None

    ###
    def __init__(self, fn=None):

        self.resdic = {}
        self.sequence = []
        self.contacts = {}
        self.stackings = {}

        if fn: self.read(fn)

        # print self.resdic
        # print self.sequence
        # print self.contacts
        # print self.stackings
        
        return None

    ###
    def save(self, fn):
        store_d = {'seq': self.sequence, \
                   'contacts': self.contacts, \
                   'stackings': self.stackings}
        fo = open(fn + '.dat', 'wb')
        pickle.dump(store_d, fo, 1)
        fo.close()
	return None


    pass # class-pass :P

###
def compute_all():
	

    all_mca = {}
    if os.path.exists('all_mca.dat'):
	    all_mca = pickle.load(open('all_mca.dat'))
    for fi in glob.glob('*.mca'):
        print fi
        mcad = MCA_Data(fi)
        all_mca[fi[:fi.rfind('.mca')]] = {'seq': mcad.sequence, \
                                          'contacts': mcad.contacts, \
                                          'stackings': mcad.stackings}
        pass
    fo = open('all_mca.dat', 'wb')
    pickle.dump(all_mca, fo, 1)
    fo.close()
        
        


###
def main(argv):

    #mcad = MCA_Data(argv[0])
    #mcad.save(argv[0])
    compute_all()

    return 0

if __name__ == '__main__': main(sys.argv[1:])
