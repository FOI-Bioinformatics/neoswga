import multiprocessing
import pickle
import os
import warnings
import numpy as np
import pandas as pd
from collections import Counter
import random
import math

cpus = int(multiprocessing.cpu_count())
min_fg_freq=float(1/100000)
max_bg_freq=float(1/150000)
min_tm=15
max_tm=45
max_gini=0.6
max_primer=500
min_amp_pred = 5
max_dimer_bp = 3
max_self_dimer_bp = 4
mismatch_penalty = 4
data_dir=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data/project')

char_to_int_dict = {'A':0, 'G':1, 'T':2, 'C':3}

def create_pool(func, input_list, cpus, initializer=None, initargs=()):
    """Calculate the results of func applied to each value in input_list using a multiprocessing pool.

    Uses context manager for proper resource cleanup (prevents zombie processes).

    Args:
        func (func): Function to apply to each value in the input list.
        input_list (list): List of values to which func will be applied.
        cpus (int): The number of processes to use.
        initializer (func): Optional function to run once per worker at startup.
        initargs (tuple): Arguments to pass to the initializer function.

    Returns:
        results: List of the results of function applied to the inputlist.
    """
    # Use context manager for automatic cleanup (close + join)
    with multiprocessing.Pool(processes=cpus, initializer=initializer, initargs=initargs) as pool:
        results = pool.map(func, input_list)
    return results


def create_pool_with_progress(func, input_list, cpus, desc="Processing", initializer=None, initargs=()):
    """Calculate results with tqdm progress bar.

    Args:
        func (func): Function to apply to each value in the input list.
        input_list (list): List of values to which func will be applied.
        cpus (int): The number of processes to use.
        desc (str): Description for progress bar.
        initializer (func): Optional function to run once per worker at startup.
        initargs (tuple): Arguments to pass to the initializer function.

    Returns:
        results: List of the results of function applied to the inputlist.
    """
    try:
        from tqdm import tqdm
        with multiprocessing.Pool(processes=cpus, initializer=initializer, initargs=initargs) as pool:
            results = list(tqdm(pool.imap(func, input_list), total=len(input_list), desc=desc))
        return results
    except ImportError:
        # Fall back to standard pool if tqdm not available
        return create_pool(func, input_list, cpus, initializer, initargs)

def flatten(l):
    """Flattens a multidimensional array.

    Args:
        l: The list to flatten.

    Returns:
        flattened_list: The flattened array.
    """

    return [item for sublist in l for item in sublist]

def mergeArrays(arr1, arr2):
    """Merges two sorted arrays, maintaining sorted order.

        Args:
            arr1: First array in sorted order.
            arr2: Second array in sorted order.

        Returns:
            arr3: The merged array in sorted order.
    """
    n1 = len(arr1)
    n2 = len(arr2)

    arr3 = []
    i = 0
    j = 0

    # Traverse both array
    while i < n1 and j < n2:
        if arr1[i] < arr2[j]:
            arr3.append(arr1[i])
            i = i + 1
        elif arr1[i] > arr2[j]:
            arr3.append(arr2[j])
            j = j + 1
        else:
            j = j + 1

    # Store remaining elements of first array
    while i < n1:
        arr3.append(arr1[i]);
        i = i + 1

    # Store remaining elements of second array
    while j < n2:
        arr3.append(arr2[j]);
        j = j + 1
    return arr3

def softmax(x):
    """Compute softmax values for each sets of scores in x.

        Args:
            x: Input value for softmax.

        Returns:
            s: Value of softmax function of x.
    """
    e_x = np.exp(x)
    return e_x / e_x.sum()

def sigmoid(x):
    """Compute the standard logistic sigmoid, mapping a real number to (0, 1).

    Args:
        x: Input value.

    Returns:
        The sigmoid of x, equal to 1 / (1 + exp(-x)).
    """
    return 1 / (1 + math.exp(-x))

def output_to_df(df, sheet_name, xls_path):
    """Output a pandas dataframe to an Excel spreadsheet.

    Uses modern pandas ExcelWriter context manager pattern.

    Args:
        df: The pandas dataframe to be written to the spreadsheet.
        sheet_name: The sheet name of the excel file.
        xls_path: The path to the excel file.
    """
    with pd.ExcelWriter(xls_path, engine='openpyxl', mode='a',
                        if_sheet_exists='replace') as writer:
        df.to_excel(writer, sheet_name=sheet_name)

def get_seq_length(genome):
    """
    Computes the sequence length of a genome in one fasta file.

    Args:
        genome: Fasta file of the genome.

    Returns:
        l: Total length of characters in the fasta file.
    """
    return sum(1 for _ in read_fasta_file(genome))

def get_all_seq_lengths(fname_genomes=None, cpus=8):
    """
    Computes the sequence length of multiple genomes.

    Args:
        fname_genomes: Fasta file of the genome.
        cpus: Number of cpus to be used.

    Returns:
        l: List of lengths of the individual genomes in all the fasta files.
    """
    # seq_lengths = create_pool(get_seq_length, fname_genomes, cpus)
    seq_lengths = []

    for fname_genome in fname_genomes:
        seq_lengths.append(get_seq_length(fname_genome))
    return seq_lengths

def longest_char_repeat(s, char):
    """Find the longest consecutive run of a given character in a string.

    Args:
        s: Input string to search.
        char: The character whose longest run is sought.

    Returns:
        Length of the longest consecutive run of ``char`` in ``s``.
    """
    max_count = 0
    if s[0] == char:
        count = 1
    else:
        count = 0

    for i,c in enumerate(s):
        if i == 0:
            continue
        if c == char and s[i-1] == char:
            count += 1
        elif c == char:
            count = 1
        if count > max_count:
            max_count = count
    return max(count, max_count)

complement_dic ={'A':'T','T':'A','G':'C','C':'G', " ": " "}
_COMPLEMENT_TABLE = str.maketrans('ATGCatgc', 'TACGtacg')
def complement(text):
    """Return the Watson-Crick complement of a DNA sequence (A<->T, G<->C).

    Case is preserved. Non-ATGC characters are left unchanged.

    Args:
        text: DNA sequence string.

    Returns:
        Complementary sequence with the same length and case as the input.
    """
    return text.translate(_COMPLEMENT_TABLE)

def get_num_mismatches(x, y):
    """Count positional mismatches between two equal-length sequences.

    Args:
        x: First sequence.
        y: Second sequence (must be the same length as ``x``).

    Returns:
        Number of positions where the two sequences differ.
    """
    num_mismatches = 0
    for i in range(len(x)):
        if x[i] != y[i]:
            num_mismatches += 1
    return num_mismatches

def longest_common_substring(s1, s2):
    """Return the length of the longest common substring.

    Uses a single-row DP approach for O(min(n,m)) space instead of O(n*m).
    """
    # Ensure s2 is the shorter string for minimal memory usage
    if len(s1) < len(s2):
        s1, s2 = s2, s1
    n, m = len(s1), len(s2)
    prev = [0] * (m + 1)
    longest = 0
    for i in range(1, n + 1):
        curr = [0] * (m + 1)
        for j in range(1, m + 1):
            if s1[i - 1] == s2[j - 1]:
                curr[j] = prev[j - 1] + 1
                if curr[j] > longest:
                    longest = curr[j]
        prev = curr
    return longest

def read_fasta_file(fname):
    """
    Read FASTA file character by character (supports .fasta, .fasta.gz).

    Uses genome_io module for automatic compression detection.
    """
    import neoswga.core.genome_io as genome_io

    # Use genome_io for automatic gzip detection
    loader = genome_io.GenomeLoader()
    for sequence in loader.load_genome_streaming(fname):
        for ch in sequence:
            yield ch

def reverse(seq):
    """Reverse a string.

    For the reverse complement of a DNA sequence, use
    :func:`reverse_complement` instead.

    Args:
        seq: Input string.

    Returns:
        The input string in reversed order.
    """
    return seq[::-1]

def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.

    DEPRECATED: Use neoswga.core.thermodynamics.reverse_complement instead.
    This function is kept for backwards compatibility.

    Args:
        seq: DNA sequence string

    Returns:
        Reverse complement of the sequence
    """
    # Use the robust implementation from thermodynamics
    from neoswga.core.thermodynamics import reverse_complement as thermo_rc
    return thermo_rc(seq)

def intersection(lst1, lst2):
    """Return elements from ``lst1`` that also appear in ``lst2``.

    Order of elements in ``lst1`` is preserved. Duplicates in ``lst1`` are
    retained if they appear in ``lst2``.

    Args:
        lst1: Primary list whose order is preserved.
        lst2: Reference list used for membership testing.

    Returns:
        List of elements from ``lst1`` present in ``lst2``.
    """
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def gini_exact(array):
    """Calculate the Gini coefficient of a numpy array. Based on bottom equation from
    http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    All values are treated equally, arrays must be 1d

    Args:
        array (array): One-dimensional array or list of floats or ints.

    Returns:
        gini_index: The gini index of the array given.
    """
    if len(array) == 0:
        return 0
    array = np.array(array, dtype=float)
    if np.amin(array) < 0:
        # Values cannot be negative:
        array -= np.amin(array)
    # Values cannot be 0:
    array += 0.0000001
    # In-place sort avoids allocating a copy
    array.sort()
    # Index per array element:
    index = np.arange(1, array.shape[0] + 1)
    # Number of array elements:
    n = array.shape[0]
    # Gini coefficient:
    return ((np.sum((2 * index - n - 1) * array)) / (n * np.sum(array)))

def most_frequent(list):
    """Return the mode (most frequent element) of a list.

    When multiple elements share the highest frequency, one is selected
    at random.

    Args:
        list: Input list of hashable elements.

    Returns:
        A single element that occurs most frequently in the input.
    """
    occurence_count = Counter(list)
    high_freq_primer = occurence_count.most_common(1)[0][0]

    high_freq_primers = []

    for k,v in occurence_count.items():
        if v == occurence_count[high_freq_primer]:
            high_freq_primers.append(k)

    return random.choice(high_freq_primers)


def get_positional_gap_lengths(positions, circular=True, seq_length=None):
    """
    Calculate gaps between consecutive positions.

    Args:
        positions: Array or list of positions (sorted or unsorted)
        circular: Whether genome is circular (wraps around)
        seq_length: Genome length (required if circular=True)

    Returns:
        Array of gap lengths between consecutive positions
    """
    if positions is None or len(positions) == 0:
        return np.array([])

    positions = np.array(positions)
    if len(positions) == 1:
        return np.array([])

    # Sort positions
    positions = np.sort(positions)

    # Calculate gaps between consecutive positions
    gaps = np.diff(positions)

    # Handle circular genome - add gap from last to first position
    if circular and seq_length is not None and seq_length > 0:
        wrap_gap = seq_length - positions[-1] + positions[0]
        gaps = np.append(gaps, wrap_gap)

    return gaps