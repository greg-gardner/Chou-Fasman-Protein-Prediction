# aminotables.py -- tables for amino acids for Chou Fasman algorithm (aka CFM)

############################# ORIGINAL PROPENSITIES #########################3
# Original CFM amino acid conformational parameters
# from [Biochemistry, Vol.13 No.2, 1974] copywright ACS
alpha_param= {  # parameters of alpha helix for CFM (Chou FasMan)
    'A':1.45,
    'R':0.79,
    'N':0.73,
    'D':0.98,
    'C':0.77,
    'Q':1.17,
    'E':1.53,
    'G':0.53,
    'H':1.24,
    'I':1.00,
    'L':1.34,
    'K':1.07,
    'M':1.2,
    'F':1.12,
    'P':0.59,
    'S':0.79,
    'T':0.82,
    'W':1.14,
    'Y':0.61,
    'V':1.14,
    'B':0.86, # average of Aspartate and Asparagine
    'Z':1.35, # average of Glutamine and Glutamate
    'X':1.00  # average for all amino acids
}
alpha_inner_param= {  # parameters for inner helix (meaning three residues 
    'A':1.59,     # on the N- and C-terminals are omitted).
    'R':0.67,
    'N':0.53,
    'D':0.53,
    'C':0.33,
    'Q':0.98,
    'E':1.45,
    'G':0.53,
    'H':0.87,
    'I':1.22,
    'L':1.91,
    'K':1.13,
    'M':1.25,
    'F':1.14,
    'P':0.0,
    'S':0.7,
    'T':0.75,
    'W':1.33,
    'Y':0.58,
    'V':1.42,
    'B':0.53, # average of Aspartate and Asparagine
    'Z':1.21, # average of Glutamine and Glutamate
    'X':0.95  # average for all amino acids
}
beta_param= {  # parameters for beta regions
    'A':0.97,
    'R':0.9,
    'N':0.65,
    'D':0.8,
    'C':1.3,
    'Q':1.23,
    'E':0.26,
    'G':0.81,
    'H':0.71,
    'I':1.6,
    'L':1.22,
    'K':0.74,
    'M':1.67,
    'F':1.28,
    'P':0.62,
    'S':0.72,
    'T':1.2,
    'W':1.19,
    'Y':1.29,
    'V':1.65,
    'B':0.73, # average of Aspartate and Asparagine
    'Z':0.75, # average of Glutamine and Glutamate
    'X':1.04  # average for all amino acids
}
coil_param= {  # parameters for coil regions
    'A':0.66,
    'R':1.2,
    'N':1.33,
    'D':1.09,
    'C':1.07,
    'Q':0.79,
    'E':0.87,
    'G':1.42,
    'H':0.92,
    'I':0.78,
    'L':0.66,
    'K':1.05,
    'M':0.61,
    'F':0.81,
    'P':1.45,
    'S':1.27,
    'T':1.05,
    'W':0.82,
    'Y':1.19,
    'V':0.66,
    'B':1.21, # average of Aspartate and Asparagine
    'Z':0.83, # average of Glutamine and Glutamate
    'X':0.99  # average for all amino acids
}
########################## HYDROPHOBICITY ###################################
# Hydrophobicity values taken from:
# Mandell AJ, Selz KA, Shlesinger MF.;
# "Wavelet transformation of protein hydrophobicity sequences suggests their memberships in structural families"
# Physica A. 1997 p.254-262
hydrophobicity= {  
    'A':0.87,
    'R':0.85,
    'N':0.09,
    'D':0.66,
    'C':1.52,
    'Q':0,
    'E':0.67,
    'G':0,
    'H':0.87,
    'I':3.15,
    'L':2.17,
    'K':1.64,
    'M':1.67,
    'F':2.87,
    'P':2.77,
    'S':0.07,
    'T':0.07,
    'W':3.77,
    'Y':2.76,
    'V':1.87,
    'B':0.38,    # average of N and D
    'Z':0.34,    # average of Q and E
    'X':1.42     # average for all amino acids
}

############################# OTHER ##########################################
# Dictionary of IUPAC amino acid abreviations per NCBI.
aminoAbrevs = {
    'Ala':'A', 'Alanine':'A',
    'Arg':'R', 'Arginine':'R',
    'Asn':'N', 'Asparagine':'N',
    'Asp':'D', 'Aspartate':'D',
    'Cys':'C', 'Cysteine':'C',
    'Gln':'Q', 'Glutamine':'Q',
    'Glu':'E', 'Glutamate':'E',
    'Gly':'G', 'Glycine':'G',
    'His':'H', 'Histidine':'H',
    'Ile':'I', 'Isoleucine':'I',
    'Leu':'L', 'Leucine':'L',
    'Lys':'K', 'Lysine':'K',
    'Met':'M', 'Metionine':'M',
    'Phe':'F', 'Phenylalanine':'F',
    'Pro':'P', 'Proline': 'P',
    'Ser':'S', 'Serine':'S',
    'Thr':'T', 'Threonine':'T',
    'Trp':'W', 'Tryptophan':'W',
    'Tyr':'Y', 'Tyrosine':'Y',
    'Val':'V', 'Valine':'V',
    'Asx':'B', # Aspartate or Asparagine
    'Glx':'Z', # Glutamine or Glutamate
    'Xaa':'X'  # Any amino acid
    }

def average_params(params):
    # -- Average up to 20 of the residues.
    if (len(params) > 20):
        print "average_params():  ERROR! Param dictionary contains more than 20 residues."
    keys = params.keys()
    total = 0
    for k in keys:
        total += float(params[k])
    print total/20

# acids = {
#     'A':,
#     'R':,
#     'N':,
#     'D':,
#     'C':,
#     'Q':,
#     'E':,
#     'G':,
#     'H':,
#     'I':,
#     'L':,
#     'K':,
#     'M':,
#     'F':,
#     'P':,
#     'S':,
#     'T':,
#     'W':,
#     'Y':,
#     'V':,
#     'B':,    # average of N and D
#     'Z':,    # average of Q and E
#     'X':     # average for all amino acids
#}