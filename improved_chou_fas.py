#!/usr/bin/python

# chou_fas.py  
# -- by  Greg Gardner  10/16/12
# -- Guesses the secondary structure of a sequence of amino acids (aka residues)
# -- using statisticaly derived conformation parameters.
# -- Based on the original 1974 paper by Peter Chou and Gerald Fasman published
# -- in [Biochemistry, Vol.13 No.2, 1974] copywright ACS

# -- Requirements:
# --   + Any shell that supports ANSII color escape sequences. 
# --     Some shells *cough Windows* require you to toggle color support.
# --   + Should work in Python 2.6 (I tested with 2.7.3).

#! In is_nucleation_region(), the minimum parameter values for determining conformations
#! should be double checked.

#! Missing predictions for 'alpha_inner' and 'coil' conformations.
 
#! In extend_regions(), does the window extend from the region or shift?
#! EX:  For sequence : [------|-|------] do you average
#!                     [----|---|------] or
#!                     [----|-|--------] ?

#! Beta sheets way over predicted.

import sys
import aminoTables
import time         #!! for testing

def main():
    SEQUENCE = testSeq[0]
#   testing()
    #generate list of indexes of regions that possess 4 of 6 residues p-Alpha() > 1.03
#    alpha_regions = find_nucleation_regions( SEQUENCE, 'alpha', 4, 6, 1.03)
#    print_regions(SEQUENCE,alpha_regions,'red')#!
#    alpha_regions_extended = extend_regions(SEQUENCE, alpha_regions, 'alpha')
#    print_regions( SEQUENCE,  alpha_regions_extended, 'blue')
    #generate list beta-regions[] of region[4] that possess 3 of 5 residues p-Beta() > 1.05  
    beta_regions = find_nucleation_regions( SEQUENCE, 'beta', 3, 5, 1.5)
    beta_regions_extended = extend_regions( SEQUENCE, beta_regions, 'beta')
    print_regions( SEQUENCE, beta_regions_extended, 'blue')

#! ############################ TESTING ####################################
def testing():#!
    return false
TIME = 0.1 #!
def wait():#!
    time.sleep(TIME)
testSeq = [
    # Heavy alpha-helix protein: Bovine RHODOPSIN (PDB# 1GZM)
    'XMNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA'
    ]

################### FUNCTIONS: CONFORMATION REGIONS ########################

def find_nucleation_regions (sequence, conformation, window_threshold, window_size, min_param):
    # -- Return list of non-redundant indexes representing nucleation
    # -- in 'sequence'. Each region contains at least 'window_threshold' residues
    # -- with 'conformation' parameters > 'min_param'.
    nucleation_regions = []
    temp_residues = []
    i = 0
    j = window_size
    while j < len(sequence):
    # Examine each window of 'window_size' in 'sequence'
        wait()#!
        window = sequence[i:j]
        if is_nucleation_region( window, conformation, window_threshold, min_param):
            nucleation_regions.append( (i,j) )
            print_regions(sequence,[(i,j)],'red')
        else: print_regions(sequence,[(i,j)],'blue')
        # shift window indexes
        i += 1
        j += 1
    #!return merge_overlapping_indexes( nucleation_regions)
    return nucleation_regions#!

def extend_regions(seq, regions, conformation):
    # -- Return a list of non-overlapping indexes representing each
    # -- of the extended 'regions' of 'seq'.
    extended_regions = []
    for reg in regions:        
    # Extend each 'reg' in 'regions':
        N = 0 # track number of shifts towards 'N' terminus
        C = 0 # track number of shifts towards 'C' terminus
        shift = (reg[0],reg[0]+4)  # slice of first four residues in 'reg'
        while(True):
        # Examine region of four shifted towards the N terminus:
            if (shift
                and average_param( conformation, seq[shift[0]:shift[1]]) >= 1 ):
                N += 1
                shift = shift_window(seq, shift, 'N') # shift towards N-term
            else:
                break
        shift = (reg[1]-4,reg[1]) # last four residues in 'reg'
        while(True):
        # Examine region of four shifted towards the C terminus:
            if (shift 
                and average_param( conformation, seq[shift[0]:shift[1]]) >= 1 ):
                C += 1
                shift = shift_window(seq, shift, 'C') # shifted towards C-term
                #!print_regions(seq, [shift], 'red')
                #time.sleep(0.5)#!
            else:
                break
        extended_regions.append((reg[0]-N, reg[1]+C)) 
    print merge_overlapping_indexes( extended_regions)
    #!return merge_overlapping_indexes( extended_regions)
    return extended_regions #!

############## FUNCTIONS:  CONFORMATIONAL PARAMETERS ######################

def is_nucleation_region( window, conformation, window_threshold, min_param):
    # -- Returns True IF 'window' contains 'window_threshold' or more of 
    # -- residues with a 'conformation' parameter > 'min_param',
    # -- AND if the average parameter value of 'window' is >= the average
    # -- parameter values for other conformation types.
    qualified = 0    
    for res in window:        
        if conformation_parameter(conformation, res) > min_param:
            qualified += 1
    if (qualified >= window_threshold
        and conformation_has_higest_average(conformation, window)):
        return True
    else: return False

def conformation_has_higest_average(conformation, sequence):
    # -- Return True if the average the average parameters of
    # -- 'conformation' in 'window' are higher than the averages
    # -- for the other three conformation types.
    average         = average_param(conformation, sequence)
    alpha_average   = average_param('alpha', sequence)
    a_inner_average = average_param('alpha_inner', sequence)
    beta_average    = average_param('beta', sequence)
    coil_average    = average_param('coil', sequence)
    has_highest_average = {
    'alpha':(            
        average > a_inner_average and
        average > beta_average and
        average > coil_average ),
    'alpha_inner':(            
        average > alpha_average and
        average > beta_average and
        average > coil_average),
    'beta':(            
        average > alpha_average and
        average > a_inner_average and
        average > coil_average ),
    'coil':(            
        average > alpha_average and
        average > a_inner_average and
        average > beta_average )
    }
    return has_highest_average[conformation]

def average_param (conformation, window):
    total = 0
    for res in window:
        total += conformation_parameter(conformation, res)
    return total/len(window)

def conformation_parameter (conformation, res):
    case = {
        'alpha':       aminoTables.alpha_param[res],
        'alpha_inner': aminoTables.alpha_inner_param[res],
        'beta':        aminoTables.beta_param[res],
        'coil':        aminoTables.coil_param[res]
        }
    return case[conformation]

################ Printing and Highlighting Sequences ########################
def print_regions(seq,regions,color):
    # -- Print 'seq' with ASCII colored 'regions'.
    higlighted_seq = ''
    if (regions):
        higlighted_seq = print_higlighted_regions(seq,regions,color)
    else: print "print_sequence(): region empty"

class colors:
    # -- ASCII color escape codes for OS shells.
    HEADER = '\033[95m'
    BLUE =   '\033[94m'
    GREEN =  '\033[92m'
    WARNING ='\033[93m'
    FAIL =   '\033[91m'
    ENDC =   '\033[0m'
    len_color = 8
    len_end =   7    

def print_higlighted_regions(seq,regions,color):
    # -- Inserts  ASCII color escape codes into 'seq'. 
    # -- Start and End codes wrap the regions of 'seq'
    # -- specified by each index pair in 'regions'.
    #! OPTIMIZE: r = index_in_regions(i, regions)
    if not regions: 
        print seq ; return
    higlighted_seq = ''
    i = 0
    len_colored = 0 
    while(i < len(seq)):
        r = index_in_regions(i,regions)
        if r:
            if (r[0] == r[1]):
                print ERRORS[1]; return           
            reg = seq[r[0]:r[1]]
            colored_reg = string_color_wrap( reg,color)
            higlighted_seq += colored_reg
            i += len(reg)
            #! need to remove region r later
        else:
            higlighted_seq += seq[i]
            i += 1
            len_colored = 0
    print higlighted_seq
    return higlighted_seq

def string_color_wrap(string, color):
    if len(string) > 0:
        if color == "error":
            return colors.FAIL + string + colors.ENDC            
        if color == "blue":
            return colors.BLUE + string + colors.ENDC
        if color == "green":
            return colors.GREEN + string + colors.ENDC
        if color == "red":
            return colors.FAIL + string + colors.ENDC

########################## Index Management ###############################
def index_in_regions(ind, regions):
    for reg in regions:
        if (ind >= reg[0] and ind < reg[1]):
            return reg
    return False        

def shift_window(seq, window, terminus):
    # -- Return the 'window' indexes in 'seq' after it is 
    # -- shifted in the direction of 'terminus'.
    if (terminus == 'N'
        and window[0]-1 >= 0):
        return (window[0]-1, window[1]-1)
    if (terminus == 'C'
        and window[1]+1 < len(seq)):
        return (window[0]+1, window[1]+1)
    return False

def merge_overlapping_indexes(indexes): 
    # -- Return list of non-overlapping index pairs from a sorted list.
    # -- Yeah that's right: 'indexes' not 'indeces'. To hell with proper grammar.
    merged_indexes = []
    if len(indexes) < 2:
        return indexes
    i = 0
    while i < len(indexes)-1:
    # step through indexes
        current = indexes[i]
        next =    indexes[i+1]
        if (merged_indexes
            and try_merging( merged_indexes[-1], current)):
            # If the last merged pair and the current pair overlap:
            merged_indexes[-1] = try_merging( merged_indexes[-1], current)
        elif try_merging( current, next):
            # If the current and next index pairs overlap:
            merged_indexes.append( try_merging( current, next))
            i += 1
        else:
            merged_indexes.append(indexes[i])
        i += 1
    return merged_indexes
def try_merging(i,j):
    # -- NOTE: these index pairs represent a Python slice!
    if (i[1] >= j[0]):
        return (i[0],j[1])
    else: return False

####################### Error Messages (colored) ###########################
ERRORS = []
for Error in [
 "print_higlighted_regions:  ERROR 0: sequence empty OR no parts of sequence are conformation regions.\n",#! not used
 "print_higlighted_regions:  ERROR 1: 'regions' contains a redundant index pair. Aborting loop.",
 ""
 ]: ERRORS.append(string_color_wrap(Error,"error"))


# Call main()
if __name__ == '__main__':
    main()
