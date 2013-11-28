#!/usr/bin/python

# chou_fas.py  
# -- by  Greg Gardner  10/16/12
# -- Guesses the secondary structure of a sequence of amino acids (aka residues)
# -- using statisticaly derived conformation parameters.
# -- Based on the original 1974 paper by Peter Chou and Gerald Fasman published
# -- in [Biochemistry, Vol.13 No.2, 1974] copywright ACS

# -- Requirements:
# --  + Any shell that supports ANSII color escape sequences. 
# --    Some shells (Windows) require you to toggle color support.
# --  + Should work in Python 2.6 (I tested with 2.7.3).


# -- Issues:
# --  + In is_nucleation_region(), the minimum parameter values for determining conformations
# --    should be double checked.
# --  + Missing predictions for 'alpha_inner' and 'coil' conformations.
# --  + In extend_regions(), does the window extend from the region or shift?
# --    EX:  For sequence : [------|-|------] do you average
# --                        [----|---|------] or
# --                        [----|-|--------] ?
# --  + Beta sheets way over predicted.

# Libraries
import sys
import amino_tables
import time         #!! for testing

# Global SEQUENCE
SEQENCE =[]

def main():
    global SEQUENCE
    
    # Print initial message.
    SEQUENCE = testSeq[0]
    user_input = 0
    print highlighted_string("chou-fas.py \nThe original Chou-Fasman Protein Conformation Prediction Algorithm", "green");
    print "\nSequence: ";  print highlighted_string(SEQUENCE,"green"); print "\nProcessing...\n"

    # Calculate the possible nucleation regions for the
    # 3 conformation types: 'alpha', 'beta', 'turn'
    alpha_regions = find_nucleation_regions ( 4, 6, 1, "alpha")
    beta_regions =  find_nucleation_regions ( 3, 5, 1, "beta")
    turn_regions =  [] #!
    
    # Extend each proposed nucleation  region in the 3' and 5' directions
    # until the average parameter value < 1. 
    alpha_regions = extend_regions( alpha_regions, "alpha")
    beta_regions =  extend_regions( beta_regions,  "beta")

    # Filter out the extended regions with average propensity values that
    # are lower than the given threshholds.
    alpha_regions = filter_extended_regions( alpha_regions, "alpha", 1.03)
    beta_regions  = filter_extended_regions( beta_regions,  "beta",  1.05)
    #turn_regions  = filter_extended_regions(turn_regions,  "turn",  1.00)

    # Merge any the extended regions that overlap. 
    alpha_regions = merge_overlapping_regions( alpha_regions)
    beta_regions =  merge_overlapping_regions( beta_regions)    
    
    # When regions of different conformation types overlap, keep the
    # region with the highest average propensity.
    alpha_regions_final = solve_overlaps( alpha_regions, 'alpha',\
                                    beta_regions, 'beta')
   # beta_regions  = solve_overlaps( beta_regions, turn_regions)
    beta_regions_final  = solve_overlaps( beta_regions, 'beta',  \
                                    alpha_regions, 'alpha')
    #turn_regions  = solve_overlaps(turn_regions, alpha_regions)

    #! Merge all regions into one highlighted sequence
#    all_regions = highlight_all_conformations(alpha_regions, beta_regions, turn_regions)
    #..............
    
    # Prompt user loop.
    while user_input >= 0:
        user_input = input(prompt_string)
        if user_input == 0:
            return
        if user_input == 1:
	    print "Alpha regions:"
            print highlighted_regions( alpha_regions_final, 'red'); 
            print "Beta regions:"
	    print highlighted_regions( beta_regions_final, 'blue'); print
        if user_input == 2:
            animate_regions(SEQUENCE, alpha_regions, "red", .05)
        if user_input == 3:
            animate_extend_regions(SEQUENCE, alpha_regions_final,"alpha", .05)
#! ############################ TESTING ####################################
TIME = .1 #!
def wait():#!
    time.sleep(TIME)
testSeq = [
    # Heavy alpha-helix protein: Bovine RHODOPSIN (PDB# 1GZM)
    'XMNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA'
    ]
prompt_string = 'Type a number and press enter:\n 1 : Display the final output of the program. \n 2 : Animate find_nucleation_regions() (Incomplete)\n 3 : Animate extend_regions() (Incomplete)\n'

################### FUNCTIONS: CONFORMATION REGIONS ########################

def find_nucleation_regions (window_threshold, window_size, min_param, conformation):
    # -- Return list of indexes representing the probable nucleation regions
    # -- in SEQUENCE. Each region contains at least 'window_threshold' residues
    # -- with 'conformation' parameters > 'min_param'.
    # -- Some regions may overlap.
    global SEQUENCE
    nucleation_regions = []
    temp_residues = []
    i = 0
    j = window_size
    # Examine each window of 'window_size' in 'SEQUENCE'.
    # If region is found, append it to 'nucleation_regions'.
    while j <= len(SEQUENCE):
        if is_nucleation_region( (i,j), conformation, window_threshold, min_param):
            nucleation_regions.append( (i,j) )
        # shift window indexes
        i += 1
        j += 1
    return nucleation_regions

def extend_regions(regions, conformation):
    # -- Return a list of indexes representing each
    # -- of the extended 'regions' in the global SEQUENCE.
    # -- Some regions may overlap.
    extended_regions = []
    for reg in regions:
    # Extend each 'reg' in 'regions':
        N = 0 # track number of shifts towards 'N' terminus
        C = 0 # track number of shifts towards 'C' terminus
        shift = (reg[0],reg[0]+4)  # slice of first four residues in 'reg'
        while(True):
        # Examine region of four shifted towards the N terminus:
            if (shift
                and average_param( shift, conformation) >= 1 ):
                N += 1
                shift = shift_window( shift, 'N') # shift towards N-term
            else:
                break
        shift = (reg[1]-4,reg[1]) # last four residues in 'reg'
        while(True):
        # Examine region of four shifted towards the C terminus:
            if (shift 
                and average_param( shift, conformation) >= 1 ):
                C += 1
                shift = shift_window( shift, 'C') # shifted towards C-term
            else:
                break
        extended_regions.append((reg[0]-N, reg[1]+C)) 
    return extended_regions

############## FUNCTIONS:  CONFORMATIONAL PARAMETERS ######################

def is_nucleation_region( region, conformation, window_threshold, min_param):
    # -- Returns True IF 'region' of SEQUENCE contains 'window_threshold' 
    # -- or more residues with a 'conformation' parameter > 'min_param',
    # -- AND if the average parameter value of 'region' is >= the average
    # -- parameter values for other conformation types.
    global SEQUENCE
    qualified = 0    
    window = SEQUENCE[ region[0] : region[1] ]
    for residue in window:        
        if conformation_parameter(residue, conformation) > min_param:
            qualified += 1
    if (qualified >= window_threshold
        and conformation_has_higest_average( region, conformation)):
        return True
    else: return False

def conformation_has_higest_average(region, conformation):
    # -- Return True if this 'region' of SEQUENCE has an average
    # -- 'conformation' propensity value that is higher than
    # -- the averages for the other conformation types.
    average         = average_param( region, conformation)
    alpha_average   = average_param( region, 'alpha')
    a_inner_average = average_param( region, 'alpha_inner')
    beta_average    = average_param( region, 'beta')
    coil_average    = average_param( region, 'coil')
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

def average_param (region, conformation):
    # -- Return the average 'conformation' parameter value
    # -- for the 'region' of SEQUENCE.
    # -- 'region' is a tuple of two indexes representing
    # -- a sub-sequence of SEQUENCE.
    total = 0
    sequence = SEQUENCE[ region[0] : region[1] ]
    for aminoAcid in sequence:
        total += conformation_parameter( aminoAcid, conformation)
    return total/len(sequence)

def conformation_parameter (residue, conformation):
    case = {
        'alpha':       amino_tables.alpha_param[ residue],
        'alpha_inner': amino_tables.alpha_inner_param[ residue],
        'beta':        amino_tables.beta_param[ residue],
        'coil':        amino_tables.coil_param[ residue]
        }
    return case[conformation]


################ Printing and Highlighting Sequences ########################

class colors:
    # -- ASCII color escape codes for OS shells.
    HEADER = '\033[95m'
    BLUE =   '\033[94m'
    GREEN =  '\033[92m'
    GREY =   '\033[1;30m'
    WARNING ='\033[93m'
    FAIL =   '\033[91m'
    ENDC =   '\033[0m'
    len_color = 8
    len_end =   7    

def print_regions(alpha_regions, beta_regions, turn_regions):
    # -- Print SEQUENCE with these 'regions' highlighted
    # -- separate colors.
    # -- NOTE: make sure there are no overlapps between regions.
    global SEQUENCE
    highlighted_regions = []
    i = 0
    # Rebuild the SEQUENCE string with the conformation regions
    # highlighted the appropriate colors.
    while i < len(SEQUENCE):
        a = index_in_regions( i, alpha_regions) # Returns the matching
        b = index_in_regions( i, beta_regions)  #  region or False.
        t = index_in_regions( i, turn_regions)
        if a:
            highlighted_regions.append( highlighted_region( a, 'red'))
            i = a[1]
        elif b:
            highlighted_regions.append( highlighted_region( b, 'blue'))
            i = b[1]
        elif t:
            highlighted_regions.append( highlighted_region( t, 'green'))
            i = t[1]
        else:
            i += 1
    print highlighted_regions
    return highlighted_regions
    

def highlighted_regions(regions,color):
    # -- Wrap each 'region' of 'seq' with ASCII color escape codes. 
    # -- Start and End codes wrap the regions of 'seq'
    # -- specified by each index pair in 'regions'.
    #! OPTIMIZE: r = index_in_regions(i, regions)
    global SEQUENCE
    if not regions: 
        print SEQUENCE ; return
    highlighted_seq = ''
    i = 0
    while(i < len(SEQUENCE)):
        r = index_in_regions(i,regions) #! 
        if r:
            if (r[0] == r[1]):
                print ERROR[1]; return
            reg = SEQUENCE[r[0]:r[1]]
            highlighted_seq += highlighted_string(reg, color)
            i += len(reg)
            #! need to remove region r later
        else:
            highlighted_seq += SEQUENCE[i]
            i += 1
            len_colored = 0
    return highlighted_seq

def highlighted_region(region, color):
    global SEQUENCE
    string = SEQUENCE[ region[0] : region[1]]
    return highlighted_string( string, color)

def highlighted_string(string, color):
    if len(string) > 0:
        if color == "red":
            return colors.FAIL + string + colors.ENDC
        if color == "blue":
            return colors.BLUE + string + colors.ENDC
        if color == "grey":
            return colors.GREY + string + colors.ENDC
        if color == "green":
            return colors.GREEN + string + colors.ENDC
        else:
            print ERROR[2]

#!
"""
def highlight_all_conformations(alphas, betas, turns):
    # -- Combines all regions of all types of conformation
    # -- and returns one final highlighted sequence of
    # -- amino acid residues.
    regions = merge_regions(alphas, betas)
    regions = merge_regions(regions, turns)
    final_sequence = 'no dice'
    return final_sequence
"""
########################## Index Management ###############################
def index_in_regions(ind, regions):
    for reg in regions:
        if (ind >= reg[0] and ind < reg[1]):
            return reg
    return False        

def shift_window(window, terminus):
    # -- Return the 'window' indexes in SEQUENCE after it is 
    # -- shifted in the direction of 'terminus'.
    if (terminus == 'N'
        and window[0]-1 >= 0):
        return (window[0]-1, window[1]-1)
    if (terminus == 'C'
        and window[1]+1 < len(SEQUENCE)):
        return (window[0]+1, window[1]+1)
    return False

def merge_overlapping_regions(indexes): 
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

def are_overlapping(reg1, reg2):
    #  -- Return true if regions overlap.
    result = False
    if reg1[1] >= reg2[0] or reg2[1] >= reg1[0] or \
       reg1[0] <= reg2[1] or reg2[0] <= reg1[1]:
       result = True
    return result

def filter_extended_regions(regions, conformation, threshold):
    # -- Return list of filtered 'regions' with 'conformation'
    # -- propensities >= 'threshold'.
    if not regions: return []
    filtered_regions = []
    for reg in regions:
        if average_param( reg, conformation) >= threshold:
            filtered_regions.append(reg)
    return filtered_regions

def solve_overlaps(regions1, conf1, regions2, conf2):
    # -- Returns all regions from 'regions1' that either:
    # -- 1) do not overlap with regions from 'regions2'.
    # -- 2) overlaps with a region from 'regions2',
    # --    but has a higher conformation affinity.
    solved_regions = []
    keep = ()
    for reg1 in regions1:
        for reg2 in regions2:
            if are_overlapping(reg1, reg2):
                keep = stronger_affinity(reg1, conf1, reg2, conf2)
            else:
                keep = reg1
            solved_regions.append(keep)
    return solved_regions

def stronger_affinity(reg1, conf1, reg2, conf2):
    # -- Returns the region with the strongest 'conformation'
    # -- affinity.
    # -- In case of tie 'alpha' conformations win.
    stronger = reg2
    if average_param(reg1, conf1) > average_param(reg2, conf2):
        stronger = reg1
    elif average_param(reg1, conf1) == average_param(reg2, conf2):
        if conf1 == 'alpha':
            stronger = reg1
        elif conf2 == 'alpha':
            stronger = reg2
    return stronger

def merge_all_regions_into_highlighted_sequence (alpha_reg, beta_reg, turn_reg):
    # -- A particularly ugly function that takes the
    # -- given regions (indexes into the Sequence)
    # -- and produces a final multicolored sequence
    # -- of amino acids.
    
    # Check if any list of regions is empty.
    if not alpha_reg and not regions2: return []
    if not regions1: return regions2
    if not regions2: return regions1
    # Proceed to merge.

    highlighted_sequence = []
    i1 = 0
    i2 = 0
    i3 = 0
    max1 = len(alpha_reg)
    max2 = len(beta_reg)
    max3 = len(turn_reg)
    while ia < max1 and ib and it :
        current1 = regions1[i1][0]
        current2 = regions2[i2][0]
        if current1 < current2:
            merged_regions.append(regions1[i1])
            i1 +=1
        else:
            merged_regions.append(regions2[i2])
            i2 +=1
    # Append whatever is left.
    if i1 == max1:
        merged_regions.append( regions2[i2:])
    else:
        merged_regions.append( regions1[i1:])
    return merged_regions


##########################   Animations   ###############################

def animate_regions(sequence, regions, color, seconds):
    # -- This is the same function as find_nucleation_regions() except
    # -- this animates steps in the terminal.
    # -- 'seconds' is step rate of the animation.
    if regions:
        window_size = regions[0][1] - regions[0][0]
    else:
        print "animate_regions:  No regions to animate."
        return
    window = (0, window_size)
    i = 0
    # Examine each window of 'window_size' in 'sequence'.
    # If region is found, append it to 'regions'.
    while window[1] <= len(sequence) and i < len(regions) :
        if window == regions[i]:
            print highlighted_string(sequence[window[0]:window[1]], color); print #colored_regions(sequence[window], color); print
            i += 1
        else:
            print highlighted_string(sequence[window[0]:window[1]],'grey'); print
        window = (window[0]+1, window[1]+1)  # shift window of sequence
        time.sleep(seconds)
    return merge_overlapping_regions(regions)

def animate_extend_regions(seq, regions, conformation, seconds):
    # -- This is the same function as extend_regions() except
    # -- this animates animates steps in the terminal.
    # -- 'seconds' is step rate of the animation.
    extended_regions = []
    for reg in regions:
    # Extend each 'reg' in 'regions':
        N = 0 # track number of shifts towards 'N' terminus
        C = 0 # track number of shifts towards 'C' terminus
        shift = (reg[0],reg[0]+4)  # slice of first four residues in 'reg'
        covered_region = reg
        while(True):
        # Examine region of four shifted towards the N terminus:
            if (shift
                and average_param( (shift[0],shift[1]), conformation) >= 1 ):
                N += 1
                print highlighted_string(seq[reg[0]-N:reg[1]+C],'red')
                time.sleep(seconds)
                shift = shift_window(seq, shift, 'N') # shift towards N-terminus
            else:
                break
        shift = (reg[1]-4,reg[1]) # last four residues in 'reg'
        while(True):
        # Examine region of four shifted towards the C terminus:
            if (shift 
                and average_param( (shift[0],shift[1]), conformation) >= 1 ):
                C += 1
                #print_region(seq,(reg[0]-N,reg[1]+C),'red')
                time.sleep(seconds)
                shift = shift_window(seq,shift, 'C') # shifted towards C-terminus
            else:
                break
        extended_regions.append((reg[0]-N, reg[1]+C)) 
    return merge_overlapping_indexes( extended_regions)

####################### Error Messages (colored) ###########################
ERROR = []
for Error in [
 "print_colored:  ERROR 0:  sequence empty OR no parts of sequence are conformation regions.",#! not used
 "colored_regions:  ERROR 1:  'regions' contains a redundant index pair. Aborting loop.",
 "highlighted_string: ERROR 2:  color not available!"
 ]: ERROR.append(highlighted_string(Error,"red"))


# Call main()
if __name__ == '__main__':
    main()
