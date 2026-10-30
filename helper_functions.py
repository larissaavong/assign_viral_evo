import math 

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    distance = 0
    valid_pos = 0

    for i, j in zip(reference_sequence, distant_sequence):
        if i == '-' or j == '-':
            continue

        if i != j:
            distance += 1

        valid_pos += 1

    proportion = distance / valid_pos

    return -0.75 * math.log(1 - (4/3) * proportion) * valid_pos

def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    purines = {'A', 'G'}
    pyrimidines = {'T', 'C'}

    num_transitions = 0
    num_transversions = 0
    num_valid_sites = 0

    for base1, base2 in zip(reference_sequence, distant_sequence):
        if base1 == '-' or base2 == '-' or base1 not in 'ACTG' or base2 not in 'ACTG':
            continue
        num_valid_sites += 1
        
        if base1 != base2:
            if (base1 in purines and base2 in purines) or (base1 in pyrimidines and base2 in pyrimidines):
                num_transitions += 1
            else:
                num_transversions += 1

    if num_valid_sites == 0:
        return 0.0

    prop_transition = num_transitions / num_valid_sites
    prop_transversion = num_transversions / num_valid_sites
    distance = -0.5 * math.log(1 - 2 * prop_transition - prop_transversion) - 0.25 * math.log(1 - 2 * prop_transversion) 

    return distance * num_valid_sites
