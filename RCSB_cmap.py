import os
import io
import re
import sys
import gzip
import shutil
import requests
import tempfile
import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import logging
import warnings

from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from Bio.PDB import PDBParser, Selection, MMCIFParser
from scipy.spatial import distance_matrix
from Bio.PDB.PDBExceptions import PDBConstructionWarning
logging.basicConfig(filename='cmap.log', encoding='utf-8', level=logging.WARNING)
warnings.filterwarnings('ignore', category=PDBConstructionWarning) #These will be saved in cmap.log

translate_aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

class RCSB():
    def __init__(self) -> None:
        pass
        # self.aa_translate = { "GLY":"G", "ALA":"A", "VAL":"V", "LEU":"L", "ILE":"I", "THR":"T", "SER":"S", "MET":"M", "CYS":"C", "PRO":"P",
        #                       "PHE":"F", "TYR":"Y", "TRP":"W", "HIS":"H", "LYS":"K", "ARG":"R", "ASP":"D", "GLU":"E", "ASN":"N", "GLN":"Q"}
    def Download_PDB(self, pdb_id : str, assembly_id : int, chain_id = None):

        if chain_id:
            #check for pdb format if you want to find a specific chain
            url = f"https://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId={pdb_id.upper()}&assemblyId={assembly_id}"
            url = url + f"&chainID={chain_id}"

            session = requests.Session()
            retry = Retry(connect=5, backoff_factor=0.5)
            adapter = HTTPAdapter(max_retries=retry)
            session.mount('http://', adapter)
            session.mount('https://', adapter)
            response = session.get(url)

            if response.status_code == 200:
                #Success, save the file
                print(f"{pdb_id} is Downloading...")
                return response.content, 'pdb'

        #Define the url that we are going to modify and ping for data
        #check for cif format to get the entire assembly
        url = f"https://files.rcsb.org/pub/pdb/data/assemblies/mmCIF/divided/{pdb_id[1:3].lower()}/{pdb_id.lower()}-assembly{assembly_id}.cif.gz"

        session = requests.Session()
        retry = Retry(connect=5, backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        session.mount('https://', adapter)
        response = session.get(url)

        if response.status_code == 200:
            #Success, save the file
            print(f"{pdb_id} is Downloading as cif format...")
            return response.content, 'cif'

        #both pdb and cif failed...
        return response.status_code


    def PDB_Processing(self, pdb, target, format_pdb):
        #remove anything that isn't protein
        #remove any hydrogens
        aa = list(translate_aa.keys())
        try:
            if format_pdb == 'pdb':
                structure = PDBParser().get_structure(target, pdb)
            else:
                structure = MMCIFParser().get_structure(target, pdb)
        except:
            return 'skip', None
        atom_number = 1
        aa_return = []
        # if len(structure.child_list) > 1:
        for chain in structure.get_chains():
            detach_residue_list = set()
            residues = Selection.unfold_entities(chain, 'R')
            for residue in residues:
                detach_id_list = []
                #remove residues that don't have a CA
                if 'CA' not in residue.child_dict.keys():
                    detach_residue_list.add(residue._id)
                # remove anything that is not protein
                if residue.full_id[-1][0] != ' ':
                    detach_residue_list.add(residue._id)
                if residue.resname not in aa:
                    detach_residue_list.add(residue._id)

                if residue.resname in aa:
                    aa_return.append(residue.resname)
                    for atom in residue:
                        # remove hydrogens
                        if atom.element.strip() == 'H':
                            detach_id_list.append(atom.id.strip())
                        if atom.element != 'H':
                            atom.serial_number = atom_number
                            atom_number+=1
                for atom_id in detach_id_list:
                    residue.detach_child(atom_id)
            for residue in detach_residue_list:
                chain.detach_child(residue)
        return structure, [translate_aa[key] for key in aa_return]

class CMAP():
    """
    Take in a PDB structure that was loaded using BioPython and produce a contact map
    """

    def __init__(self,structure, olig = None, chain = None, chains_like = None, leven = 30):
        #target a specific chain
        if chain:
            mono = [key for key in structure[0].child_dict.keys() if key != chain]
            for key in mono:
                structure[0].detach_child(key)

        if chains_like:
            chain_info = {}; aa_dict = {}
            for chain in structure[0].get_chains():
                residues = Selection.unfold_entities(chain, 'R')
                chain_info[chain] = sorted([(int(res._id[1]), res.resname) for res in residues], key = lambda x: x[0])
            for chain in chain_info:
                aa_dict[chain] = ''.join([translate_aa[i[1]] for i in chain_info[chain]])
                if chain.full_id[-1] == chains_like:
                    ref = ''.join([translate_aa[i[1]] for i in chain_info[chain]])
            remove = []
            for aa, chain in zip(aa_dict, structure[0].child_dict.keys()):
                score, mem = self.levenshtein_dist_w_memory(ref,aa_dict[aa])
                if score > leven:
                    remove.append(chain)
            for chain in remove:
                structure[0].detach_child(chain)
            print()
            print(f"Chain IDs that were retained: {[key for key in structure[0].child_dict.keys()]}")
            print()


        self.structure = structure
        self.Load_Data()

        #normalize homo-oligomer lengths for cmap visualization
        if olig:
            self.Check_Length()

    def Load_Data(self):
        """
        Load_Data reads in the pdb residues and fills in gaps within chains
        """
        detach_list = []
        for chain in self.structure[0].get_chains():
            residues = Selection.unfold_entities(chain, 'R')

            if len(residues) == 0:
                detach_list.append(chain)
                continue
            last_res = -1
            for residue in residues:
                if residue.get_id()[1] > last_res +1 and residue.get_id()[1] != residues[0].get_id()[1]:
                    for idx in range(last_res+1,residue.get_id()[1]):
                        new_res = residue.copy()
                        for atom in new_res:
                            atom.detach_parent()
                            atom.set_coord(np.array([9999,9999,9999], dtype="float32"))
                        new_res.id = (' ', idx, chain.get_id())
                        chain.add(new_res)
                last_res = residue.get_id()[1]
        if detach_list:
            for chain in detach_list:
                self.structure[0].detach_child(chain.id)

    def levenshtein_dist_w_memory(self, string_1, string_2):
        """
        Calculate the levenshtein distance between two strings and retain a 'memory' of the edits
        to change string_2 into string_1

        """

        len_string_1 = len(string_1) + 1
        len_string_2 = len(string_2) + 1
        leven_mtx = [[(0,'') for _ in range(len_string_2)] for _ in range(len_string_1)]
        # tuple is to keep track of the operations required to edit the string for it to match
        # initialize leven_mtx
        for col in range(len_string_1):
            leven_mtx[col][0] = (col, col*'D') #D for deletion
        for row in range(len_string_2):
            leven_mtx[0][row] = (row, row*'I') #I for insertion

        # calculate full levenshtein matrix
        for col in range(1, len_string_1):
            for row in range(1, len_string_2):
                if string_1[col -1] == string_2[row -1]:
                    cost = 0
                    edit = "M" # M for match
                else:
                    cost = 1
                    edit = "S" # S for substitution
                #store values for each possibility
                D_cost, D_edit = leven_mtx[col -1][row]
                I_cost, I_edit = leven_mtx[col][row -1]
                S_cost, S_edit = leven_mtx[col -1][row -1]

                # determine the most cost effective move and store the values
                if D_cost < I_cost and D_cost + 1 < S_cost + cost:
                    leven_mtx[col][row] = (D_cost + 1, D_edit + "D")
                elif I_cost < D_cost and I_cost + 1 < S_cost + cost:
                    leven_mtx[col][row] = (I_cost +1, I_edit + "I")
                else:
                    leven_mtx[col][row] = (S_cost +cost, S_edit + edit)

        return leven_mtx[-1][-1]

    def Check_Length(self):
        """
        Check the length of each A.A. sequence in a homo-oligomer to standardize cmap size

        """
        chain_ids = list(self.structure[0].child_dict.keys())
        if len(chain_ids) >= 2:

            def check_aa(domain):
                chain_info = {}; aa_mtx = []
                for chain in self.structure[0].get_chains():
                    residues = Selection.unfold_entities(chain, 'R')
                    chain_info[chain] = sorted([(int(res._id[1]), res.resname) for res in residues], key = lambda x: x[0])
                for chain in chain_info:
                    aa_mtx.append(''.join([translate_aa[i[1]] for i in chain_info[chain]]))

                if domain == 'NTD':
                    for ref in aa_mtx:
                        mtx = []
                        for chain,aa in zip(chain_info.keys(), aa_mtx):
                            _ = self.levenshtein_dist_w_memory(ref, aa)
                            mtx.append(_)
                            if 'I' in _[-1][:10]:
                                break
                    return mtx
                if domain == 'CTD':
                    ref = max(aa_mtx, key = len)
                    mtx = []
                    for chain,aa in zip(chain_info.keys(), aa_mtx):
                        _ = self.levenshtein_dist_w_memory(ref, aa)
                        mtx.append(_)
                    return mtx

            #check NTD
            ntd = check_aa('NTD')
            for leven, chain in zip(ntd, self.structure[0].get_chains()):
                residues = sorted(Selection.unfold_entities(chain, 'R'), key = lambda x: x._id[1])
                residue = residues[0]
                end = residues[0]._id[1]
                _ = np.array([*leven[-1]], dtype='<U1')
                beg = end - np.argmax(_ == 'M')
                for idx in range(beg,end):
                        new_res = residue.copy()
                        for atom in new_res:
                            atom.detach_parent()
                            atom.set_coord(np.array([9999,9999,9999], dtype="float32"))
                        new_res.id = (' ', idx, chain.get_id())
                        chain.add(new_res)
            #check CTD
            ctd = check_aa('CTD')
            for leven, chain in zip(ctd, self.structure[0].get_chains()):
                residues = sorted(Selection.unfold_entities(chain, 'R'), key = lambda x: x._id[1])
                residue = residues[-1]
                beg = residue._id[1] + 1
                _ = np.array([*leven[-1]], dtype='<U1')
                end = beg + np.sum(_ == 'D')
                for idx in range(beg,end):
                        new_res = residue.copy()
                        for atom in new_res:
                            atom.detach_parent()
                            atom.set_coord(np.array([9999,9999,9999], dtype="float32"))
                        new_res.id = (' ', idx, chain.get_id())
                        chain.add(new_res)

    def Create_Map(self, cutoff):
        """
        Create the contact map from the PDB structure
        gaps should be filled in and homo-oligomer chain lengths should match
        even if crystal density was not present
        """
        chains = {}
        for chain in self.structure[0].get_chains():
            residues = Selection.unfold_entities(chain, 'R')
            residues = sorted(residues, key=lambda r:r.get_id()[1])
            chains[chain] = [residue['CA'] for residue in residues if 'CA' in residue.child_dict.keys()]

        coor_sep = np.array([[atom.coord for atom in chain[1]] for chain in chains.items()],dtype=object)
        coor_idx = [[i for i,a in enumerate(chain)] for chain in coor_sep]
        tot_idx,offset = [],0
        for idx, coor in enumerate(coor_idx):
            if idx == 0:
                tot_idx.append(list(coor_idx[idx]))
            else:
                offset = offset + coor_idx[idx - 1][-1] +1
                tot_idx.append([pos+offset for pos in coor_idx[idx]])
        tot_idx_flat = [pos for chain in tot_idx for pos in chain]

        coor = np.array([pos for chain in coor_sep for pos in chain])
        dist_matrix = distance_matrix(coor[tot_idx_flat,:],coor[tot_idx_flat,:],p=2)

        #contact map
        #only get contacts within predefined cutoff excluding placeholders 
        contact_map = np.zeros((len(tot_idx_flat),len(tot_idx_flat)))
        contact_map[(dist_matrix < cutoff) & (dist_matrix != 0.)] = 1

        res_list = [res for key in chains for res in chains[key]]
        combo = list(itertools.permutations(res_list,2))
        res_idx = list(itertools.permutations(list(range(len(res_list))),2))

        #fill in where all heavy atom contacts exits
        # comment these loops out for CA contacts
        fill = []
        for key,idx in zip(combo,res_idx):
            key_a,key_b = key[0],key[1]
            if key_a.parent is not None and key_b.parent is not None and idx[1] > idx[0]:
                a = np.array([key_a.parent.child_dict[key].coord for key in key_a.parent.child_dict])
                b = np.array([key_b.parent.child_dict[key].coord for key in key_b.parent.child_dict])
                dist_temp = distance_matrix(a,b,p=2)
                if np.sum(dist_temp <= cutoff).astype(bool):
                    fill.append(idx)
        for i,j in fill:
            contact_map[i][j] = 1

        #last remove hits within +-3 of the diagonal
        mask = np.zeros((len(tot_idx_flat),len(tot_idx_flat)))
        mask = np.abs(np.arange(len(tot_idx_flat)) - np.arange(len(tot_idx_flat))[:,np.newaxis]) <= 3
        contact_map[mask] = 0

        #upper triangle is a contact map that represents every contact between any heavy atom of two residues
        #lower triangle is a contact map that represents contacts between CAs of two residues
        tu = np.triu_indices(contact_map.shape[0])
        contact_map[tu[::-1]] = contact_map[tu]

        #create a mask that makes oligomers easier to visualize
        #as is this will only handle up to 5 chains, makes interchain squares darker
        mask = np.zeros((len(tot_idx_flat),len(tot_idx_flat)))
        perm = list(range(len(tot_idx)))
        perm = [p for p in itertools.permutations(perm, r=2) if p[0] < p[1]]
        value = 0.2
        for i in perm:
            mask_idx = np.array(list(itertools.product(tot_idx[i[0]],tot_idx[i[1]])))
            mask[mask_idx[:,0],mask_idx[:,1]] = 1 - int(i[1]-i[0])*value
            mask[mask_idx[:,1],mask_idx[:,0]] = 1 - int(i[1]-i[0])*value

            #recolor any contact in these regions to separate inter from intra
            _ = mask_idx[contact_map[mask_idx[:,0],mask_idx[:,1]] == 1]
            contact_map[_[:,0], _[:,1]] = 2
            contact_map[_[:,1], _[:,0]] = 2

        return contact_map,mask

    def Collapse(self,cmap):
        """ 
        Create a superimposition of all possible contacts for a homo-oligomer
        """
        chains = {}
        chain_ids = list(self.structure[0].child_dict.keys())
        if len(chain_ids) >= 2:
            for chain in self.structure[0].get_chains():
                residues = Selection.unfold_entities(chain, 'R')
                res_len = len(residues)
                chains[chain.id] = res_len

            new_cmap = np.zeros((chains[chain_ids[0]],chains[chain_ids[0]]))

            increment=0
            idx = [0]
            for _id in chains.keys():
                idx.append(increment+chains[_id])
                increment+=chains[_id]

            intra_idx = [p for p in itertools.product(idx, repeat=2) if p[0] < p[1]]
            intra_idx = [p for p in intra_idx if p[0] == p[1]-chains[chain_ids[0]]]

            # This is now a superimposition of the intramolecular contacts not one individual instance
            intra = new_cmap
            for idx in intra_idx:
                intra = intra + cmap[idx[0]:idx[1],idx[0]:idx[1]]
            mask = np.nonzero(intra)
            intra[mask] = 1

            # Superimposition of all possible intermolecular contacts
            inter = new_cmap
            inter_idx = [p for p in itertools.product(intra_idx, repeat=2) if p[0][0]+p[0][1] > p[1][0]+p[1][1]]
            for idx in inter_idx:
                inter = inter + cmap[idx[0][0]:idx[0][1],idx[1][0]:idx[1][1]]
                inter = inter + cmap[idx[1][0]:idx[1][1],idx[0][0]:idx[0][1]]

            mask = np.nonzero(inter)
            inter[mask] = 2

            #last remove hits within +-3 of the diagonal
            mask = np.zeros(inter.shape)
            mask = np.abs(np.arange(mask.shape[0]) - np.arange(mask.shape[1])[:,np.newaxis]) <= 3
            inter[mask] = 0
            new_cmap = intra + inter
            return new_cmap

        _ = np.nonzero(cmap)
        cmap[_] = 1

        return cmap

def pad_with(vector,pad_width,iaxis,kwargs):
    """
    custom paramter to create a padded edge of a given value around an np.array using np.pad()
    """
    pad_value = kwargs.get('padder', 0)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value

def overlap_definition(A, B, mtx_return=False):
    #pad edges for sliding window consistency
    padded_A = np.pad(A, ((1,1),(1,1)), mode='constant', constant_values=0)
    #Extract all possible windows
    windows = np.lib.stride_tricks.sliding_window_view(padded_A, (3,3))
    #positions with contacts
    mask_1 = np.where(B == 1)
    mask_2 = np.where(B == 2)
    mask_3 = np.where(B == 3)
    #find number of contacts in B that match a (3,3) window in A that has a nonzero elements
    matches_1 = np.any(windows[mask_1] != 0, axis=(-1,-2))
    matches_2 = np.any(windows[mask_2] != 0, axis=(-1,-2))
    matches_3 = np.any(windows[mask_3] != 0, axis=(-1,-2))

    # return B array with +s where a +-1 match in A exists and -s where they don't exist
    # 0s where no contact exists, this function retains the type of contact
    if mtx_return == True:
        B_true = np.copy(B)
        B_true[mask_1] = matches_1 * 1
        B_true[mask_2] = matches_2 * 2
        B_true[mask_3] = matches_3 * 3

        B_false = np.copy(B)
        B_false[mask_1] = ~matches_1 * 1
        B_false[mask_2] = ~matches_2 * 2
        B_false[mask_3] = ~matches_3 * 3
        B_false = B_false*-1
        return B_true + B_false
    return sum([*matches_1, *matches_2, *matches_3])

def COMPARE(_cmap1,_cmap2,pdb1_name,pdb2_name):
    """
    Compare two contact maps and return a contact map alignment based on number of symmetric contacts
    """
    A = _cmap1
    A_name = pdb1_name
    B = _cmap2
    B_name = pdb2_name

    if B.shape[0] < A.shape[0]:
        B = np.pad(B, ((0,A.shape[0]-B.shape[0]),(0,A.shape[0]-B.shape[0])), 'constant')
       
    #pad first matrix to allow for sliding window
    n = int(B.shape[0]*0.1)
    B = np.pad(B, n, pad_with, padder=0)
    mask = np.triu(np.ones(B.shape))
    B = B * mask

    #initialize the maximum number of overlapping 1s
    max_overlap = 0

    #initialize the starting idices of the optimal alignment
    start_i,start_j = 0,0

    beg = range(B.shape[0]-A.shape[0]+1)
    end = range(B.shape[0]-A.shape[0]+1)

    for i,j in zip(beg,end):
        subset = B.copy()
        subset = subset[i:i+A.shape[0], j:j+A.shape[1]]
        # only consider exact mathches
        # overlap = len(np.argwhere(np.multiply(A, subset) == 1))

        # consider matches within +- 1 | good luck
        overlap = overlap_definition(A, subset)
        #update the maximum number of overlapping 1s and the best starting indices if needed
        if overlap > max_overlap:
            max_overlap = overlap
            start_i,start_j = i,j

    #print the maximum number of overlapping 1s and the starting indicies of the optimal alignment
    print(f"Maximum number of overlapping 1s: {max_overlap}")
    print(f"Starting indices of the optimal alignment: ({start_i}, {start_j})")

    lower_triangle_idx = np.tril_indices(A.shape[0],0)

    #offset indices to create the aligned dualfold cmap
    lower_triangle_idx_align = np.transpose(lower_triangle_idx)
    lower_triangle_idx_align = lower_triangle_idx_align + start_i
    lower_triangle_idx_align = np.transpose(lower_triangle_idx_align)
    lower_triangle_idx_align = tuple(np.array(i) for i in lower_triangle_idx_align)

    #create the duafold cmap
    B[lower_triangle_idx_align] = A[lower_triangle_idx]

    return B, B_name, n, A_name, start_i, max_overlap


def Create_DF(contact_map,y_name,x_name,name="cmap",x_offset=0,y_offset=0, mask=[], olig=False):
    """ 
    Create pandas data frame of contact map for visualizing cmap
    """

    low = np.tril(contact_map) 
    low[np.triu_indices(low.shape[0])] = low.T[np.triu_indices(low.shape[0])]

    up  = np.triu(contact_map) 
    up[np.tril_indices(up.shape[0])] = up.T[np.tril_indices(up.shape[0])]

    low = overlap_definition(up, low, mtx_return=True)
    up = overlap_definition(low, up, mtx_return=True)
    # colors = {'INTRA_COMMON':'#383838', 'INTRA_UNIQUE':'#ff63f8', 'INTER_COMMON':'#7b7b7b', 'INTER_UNIQUE':'#f9c6ff'}  #pink version
    colors = {'INTRA_COMMON':'#383838', 'INTRA_UNIQUE':'#377c2b', 'INTER_COMMON':'#7b7b7b', 'INTER_UNIQUE':'#b9dcb3'}  #green version

    #check for all contacts to find matches to ACE contacts
    contact_types = {}
    contact_types["low_true_intra"]  = np.argwhere((np.tril(low) == 1) | (np.tril(low) == 3)).tolist()
    contact_types["low_false_intra"] = np.argwhere((np.tril(low) == -1) | (np.tril(low) == -3)).tolist()
    contact_types["low_true_inter"]  = np.argwhere((np.tril(low) == 2)).tolist()
    contact_types["low_false_inter"] = np.argwhere((np.tril(low) == -2)).tolist()
    contact_types["up_true_intra"]   = np.argwhere((np.triu(up) == 1) | (np.triu(up) == 3)).tolist()
    contact_types["up_false_intra"]  = np.argwhere((np.triu(up) == -1) | (np.triu(up) == -3)).tolist()
    contact_types["up_true_inter"]   = np.argwhere((np.triu(up) == 2)).tolist()
    contact_types["up_false_inter"]  = np.argwhere((np.triu(up) == -2)).tolist()

    df = {'i': [],'j': [],'type': []}
    for key in contact_types.keys():
        for i,j in contact_types[key]:
            if key == "low_true_intra":
                df['i'].append(i)
                df['j'].append(j)
                df['type'].append("INTRA_COMMON")
            if key == "low_false_intra":
                df['i'].append(i)
                df['j'].append(j)
                df['type'].append("INTRA_UNIQUE")
            if key == "low_true_inter":
                df['i'].append(i)
                df['j'].append(j)
                df['type'].append("INTER_COMMON")
            if key == "low_false_inter":
                df['i'].append(i)
                df['j'].append(j)
                df['type'].append("INTER_UNIQUE")
            if key == "up_true_intra":
                df['i'].append(i)
                df['j'].append(j)
                df['type'].append("INTRA_COMMON")
            if key == "up_false_intra":
                df['i'].append(i)
                df['j'].append(j)
                df['type'].append("INTRA_UNIQUE")
            if key == "up_true_inter":
                df['i'].append(i)
                df['j'].append(j)
                df['type'].append("INTER_COMMON")
            if key == "up_false_inter":
                df['i'].append(i)
                df['j'].append(j)
                df['type'].append("INTER_UNIQUE")

    df = pd.DataFrame.from_dict(df)
    f, ax = plt.subplots(1,1,figsize=(9,9))
    df['i'] = df['i'] - x_offset
    df['j'] = df['j'] - y_offset
    axes = []

    #keep contacts from being on the edge of the plot
    ax.set_ylim((min(df.i.min(),df.j.min())-5,max(df.i.max(),df.j.max())+5))
    ax.set_xlim((min(df.i.min(),df.j.min())-5,max(df.i.max(),df.j.max())+5))

    #standardize size of markers no matter the size of the protein
    r = 0.5 #radius of markers
    r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0] #translate radius into image coords


    #each 'type' needs its own handle for matplotlib to give unique legend elements
    for t in df['type'].unique():
        axes.append(ax.scatter(x=df.loc[df['type'].eq(t), 'i'], y= df.loc[df['type'].eq(t), 'j'],c=df.loc[df['type'].eq(t), 'type'].map(colors), s=np.pi * r_, linewidth=0, linestyle="None"))

    if len(mask) > 0:
        axes.append(ax.imshow(mask,cmap='Greys', interpolation='none', alpha=0.11))
    ax.set_ylabel(y_name)
    ax.set_xlabel(x_name)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height* 0.85])

    # Put a legend to the right of the current axis
    ax.legend(df['type'].unique(),loc='center left', bbox_to_anchor=(1, 0.5))

    name = name.replace(re.search(r"(.cif|.pdb)", name).group(0), '') if bool(re.search(r"(.cif|.pdb)", name)) else name
    print('Writing csv...',f"{name}_df.csv")
    df.to_csv(f'{name}_df.csv')
    print('Writing image...',f"{name}.png")
    plt.savefig(f"{name}.png")

    plt.clf()
    plt.close()


def main(target, chain, olig, pdb2, chain2, chains_like, leven):

    if os.path.isfile(target):
        get_pdb = RCSB()
        structure_processed, aa_structure_processed = get_pdb.PDB_Processing(target, target[0:4], target[-3:])
    else:
        get_pdb = RCSB()
        structure, format_pdb = get_pdb.Download_PDB(target, 1, chain) #download assembly 0 from rcsb
        with tempfile.NamedTemporaryFile() as temp_file:
            if format_pdb == 'pdb':
                temp_file.write(structure)
            else:
                with gzip.open(io.BytesIO(structure), 'rb') as gz: #create an in-memory stream to read the gzipped data
                    gz_bytes = gz.read()#.decode('utf-8')          #can decode to read string but write function is expecting bytes-like-object
                temp_file.write(gz_bytes)
            # biopython pdbparser expects file_type object, tempfile is used to avoid clutter
            # remove hydrogens and non polymer atoms
            structure_processed, aa_structure_processed = get_pdb.PDB_Processing(temp_file.name, target, format_pdb)

    Cmap = CMAP(structure_processed, olig=olig, chain=chain, chains_like=chains_like, leven=leven) #cmap object now contains structural information!!!
    cmap1,mask1 = Cmap.Create_Map(cutoff=8) # cutoff in angstom


    Create_DF(cmap1,target,target,name=f"{target}", mask=mask1)
    if olig:
        if len(structure_processed[0].child_list) > 1:
            cmap1_collapse = Cmap.Collapse(cmap1)
            Create_DF(cmap1_collapse,target,target,name=f"{target}_collapse",olig=olig)
        else:
            cmap1_collapse = cmap1

    if pdb2:
        if os.path.isfile(pdb2):
            get_pdb = RCSB()
            structure_processed2, aa_structure_processed2 = get_pdb.PDB_Processing(pdb2, pdb2[0:4], pdb2[-3:])
        else:
            get_pdb = RCSB()
            structure2, format_pdb2 = get_pdb.Download_PDB(pdb2, 1, chain2) #download assembly 0 from rcsb
            with tempfile.NamedTemporaryFile() as temp_file:
                if format_pdb2 == 'pdb':
                    temp_file.write(structure2)
                else:
                    with gzip.open(io.BytesIO(structure2), 'rb') as gz: #create an in-memory stream to read the gzipped data
                        gz_bytes = gz.read()#.decode('utf-8')          #can decode to read string but write function is expecting bytes-like-object
                    temp_file.write(gz_bytes)
                # biopython pdbparser expects file_type object, tempfile is used to avoid clutter
                # remove hydrogens and non polymer atoms
                structure_processed2, aa_structure_processed2 = get_pdb.PDB_Processing(temp_file.name, pdb2, format_pdb)

        Cmap2 = CMAP(structure_processed2, olig=olig, chain=chain2) #cmap object now contains structural information!!!
        cmap2,mask2 = Cmap2.Create_Map(cutoff=8) # cutoff in angstom

        Create_DF(cmap2,pdb2,pdb2,name=f"{pdb2}", mask=mask2)
        comp,x_name,padding_offset,y_name,align_offset, max_overlap = COMPARE(cmap1, cmap2,target,pdb2)
        Create_DF(comp,x_name,y_name,name=f"comp_{target}_{pdb2}")

        if olig:
            if len(structure_processed2[0].child_list) > 1:
                cmap2_collapse = Cmap2.Collapse(cmap2)
                Create_DF(cmap2_collapse,pdb2,pdb2,name=f"{pdb2}_collapse",olig=olig)
            else:
                cmap2_collapse = cmap2
            comp,x_name,padding_offset,y_name,align_offset, max_overlap = COMPARE(cmap1_collapse, cmap2_collapse,target,pdb2)
            Create_DF(comp,x_name,y_name,name=f"comp_{target}_{pdb2}_collapse")



def CL_input():
    """
    Parse command line arguments that are being passed in
    """

    error_message = """
Missing command line arguments!
Available flags:
-pdb ####      |  RCSB pdb id for the protein of interest (example: 1fha)
               |  if ####.pdb(cif) is in cwd the local file will be used (example: 1fha.pdb)
-chain #       |  Chain id of interest (example: A)
-pdb2 ####     |  used to create a dualfold comparison contact map with -pdb
-chain2 #      |  specify the chain of -pdb_2 that will be used for the comparison
-oligomer      |  If -pdb is an oligomer this flag will create a superposition of of all contacts
-chains_like # |  retain only chains that are similar to the selected chain. Useful for hetero-oligomer (example: C)
    -leven #   |  
** NOTE: chains_like calculates the levenshtein distance between chains and retains chains that are within 30
** NOTE: of the adjust -leven if this is to restrictive/permissive
"""

    options = ['-pdb', '-chain', '-pdb2', '-chain2', '-oligomer', '-chains_like', '-leven']
    for i in sys.argv:
        if '-' in i and i not in options:
            print(f"This is not an option: {i}")
            print(error_message)
            sys.exit()

    #-pdb has to be defined
    if not any((True if _ == '-pdb' else False for _ in sys.argv)) or len(sys.argv) <= 2:
        print(error_message)
        sys.exit()

    pdb = sys.argv[[idx for idx, _ in enumerate(sys.argv) if '-pdb' == _][0] + 1]
    pdb2 = None if not [idx for idx, _ in enumerate(sys.argv) if '-pdb2' == _] else sys.argv[[idx for idx, _ in enumerate(sys.argv) if '-pdb2' == _][0] + 1]
    chain = None if not [idx for idx, _ in enumerate(sys.argv) if '-chain' == _] else sys.argv[[idx for idx, _ in enumerate(sys.argv) if '-chain' == _][0] + 1]
    chain2 = None if not [idx for idx, _ in enumerate(sys.argv) if '-chain2' == _] else sys.argv[[idx for idx, _ in enumerate(sys.argv) if '-chain2' == _][0] + 1]
    chains_like = None if not [idx for idx, _ in enumerate(sys.argv) if '-chains_like' == _] else sys.argv[[idx for idx, _ in enumerate(sys.argv) if '-chains_like' == _][0] + 1]
    leven = 30 if not [idx for idx, _ in enumerate(sys.argv) if '-leven' == _] else sys.argv[[idx for idx, _ in enumerate(sys.argv) if '-leven' == _][0] + 1]
    olig = True if any((True if _ == '-oligomer' else False for _ in sys.argv)) else False
    return pdb, pdb2, chain, chain2, olig, chains_like, leven

if __name__ == "__main__":
    """
    execute code by passing in a valid rcsb id
        python RCSB_cmap.py 1fha
    you can specify a chain id as well
        python RCSB_cmap.py 1fha A

    A collapsed version of homo-oligomers can be created by passing the flag olig
        python RCSB_cmap.py 1fha olig
    """

    pdb, pdb2, chain, chain2, olig, chains_like, leven = CL_input()
    print('-pdb =', pdb,'-pdb2 =', pdb2,'-chain =', chain,'-chain2 =', chain2,'-oligomer =', olig, '-chains_like =', chains_like, 'leven =', leven)
    main(pdb, chain, olig, pdb2, chain2, chains_like, leven)
