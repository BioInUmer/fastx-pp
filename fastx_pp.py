#!/usr/bin/python3

"""
fastx_pp.py - FASTA/FASTQ Preprocessor

A command-line tool for preprocessing FASTA and FASTQ sequence files.
Supports three main operations:
    - Reverse-complement: Generate reverse complement of sequences
    - Trim: Remove bases from both ends of sequences
    - Adaptor-removal: Remove adaptor sequences from the start of reads

Usage:
    python3 fastx_pp.py --input <file> --output <file> --operation <op> [options]

Author: Muhammad Umer Hussain Kousar
Version: 1.0
"""

# Module required for command line argument parsing
import sys
# Module for matching operations
import re

# ========== AUXILIARY FUNCTIONS ==========


def dot_number(number):
    """
    Format a number with dots as thousand separators for improved readability.

    Parameters:
        number (int or float): The numeric value to be formatted. 
            Can be any positive or negative integer or floating-point number.

    Returns:
        str: A string representation of the number with dots as thousand separators.
            Examples: 1000 → "1.000", 1000000 → "1.000.000", 150 → "150"
    """
    return f"{number:,}".replace(",", ".")


def base_percentages(base_counts, bases_processed):
    """
    Calculate the percentage composition of each nucleotide base type.

    Parameters:
        base_counts (dict): A dictionary containing counts for each base type.
            Expected keys are 'A', 'C', 'G', 'T', and 'N', with integer values
            representing the count of each base.
        bases_processed (int): The total number of bases processed. Must be greater
            than zero to avoid division by zero.

    Returns:
        tuple: A tuple of five integers (per_a, per_c, per_g, per_t, per_n)
            representing the rounded percentage of A's, C's, G's, T's, and N's
            respectively in the processed bases.
    """
    per_a = int(round((base_counts['A'] / bases_processed) * 100))
    per_c = int(round((base_counts['C'] / bases_processed) * 100))
    per_g = int(round((base_counts['G'] / bases_processed) * 100))
    per_t = int(round((base_counts['T'] / bases_processed) * 100))
    per_n = int(round((base_counts['N'] / bases_processed) * 100))
    return per_a, per_c, per_g, per_t, per_n


def usage_instructions():
    """
    Display comprehensive usage instructions and exit the program.

    This function is called when incorrect arguments are provided.
    It prints detailed usage information including all available operations
    and their required parameters, then terminates execution.

    Parameters:
        None

    Returns:
        None: This function does not return. It terminates the program via sys.exit(1).
    """
    print("\n"
          "The arguments provided are incorrect!\n"
          "Usage: python3 fastx_pp.py --input <input_file> --output <output_file> --operation <operation> [--trim-left <num>] [--trim-right <num>] [--adaptor <sequence>]\n\n"
          "Valid operations are:\n"
          " • rc                : Reverse-complement the sequences in the input FASTA/FASTQ file.\n"
          " • trim              : Trim a specified number of bases from both ends of each sequence.\n"
          "                       Requires --trim-left and --trim-right arguments.\n"
          " • adaptor-removal   : Remove a specified adaptor sequence from the beginning of each sequence.\n"
          "                       Requires --adaptor argument.\n")
    print("Example: python3 fastx_pp.py --input reads.fasta --output reads_rc.fasta --operation rc")
    sys.exit(1)


# ========== INPUT VERIFICATION ==========


def check_arguments():
    """
    This function processes command-line arguments to extract and validate:
    - Input and output file paths
    - Operation type (rc, trim, or adaptor-removal)
    - Operation-specific parameters (trim lengths, adaptor sequence)

    It performs comprehensive validation including:
    - Flag recognition and uniqueness
    - Value type checking (file extensions, integers, sequences)
    - Required parameter presence
    - File extension matching

    Parameters:
        None: Reads arguments directly from sys.argv

    Returns:
        tuple: A 6-element tuple containing:
            - input_file (str): Path to the input FASTA/FASTQ file
            - output_file (str): Path to the output FASTA/FASTQ file
            - operation (str): Operation to perform ('rc', 'trim', or 'adaptor-removal')
            - ltrim_len (int): Number of bases to trim from the left end (default: 0)
            - rtrim_len (int): Number of bases to trim from the right end (default: 0)
            - adaptor (str): Adaptor sequence to remove from reads (default: empty string)
    """

    # Initialize variables for the argumenets that will be provided on the command line
    input_file = ""
    output_file = ""
    operation = ""
    ltrim_len = 0
    rtrim_len = 0
    adaptor = ""
    seen_flags = set()  # track which flags have been seen (set() to avoid duplicates)

    valid_flags = ["--input", "--output", "--operation",
                   "--trim-left", "--trim-right", "--adaptor"]
    valid_operation_values = ["rc", "trim", "adaptor-removal"]

    # Start from index: 1 to skip the script name (index: 0)
    i = 1

    # As long as there are arguments to process
    while i < len(sys.argv):

        # Get the flag (--input, --output, --operation, --trim-left, --trim-right, --adaptor)
        flag = sys.argv[i].lower()

        # ··· FLAGS CHECK ···

        # Ensure the flag is valid
        if flag not in valid_flags:
            print(f"\nError: '{flag}' not recognized")
            print(f"Valid flags are: '{valid_flags}'")
            usage_instructions()

        # Ensure that there is a value after the last flag
        if i + 1 >= len(sys.argv):
            print(f"\nError: '{flag}' requires a value")
            usage_instructions()

        # Ensure that each flag is only provided once
        if flag in seen_flags:
            print(f"\nError: '{flag}' provided multiple times")
            usage_instructions()
        # Mark this flag as seen
        seen_flags.add(flag)

        # Get the value associated with the flag (file.fastx, rc, trim, adaptor-removal, integer, sequence)
        value = sys.argv[i + 1].strip()

        # ··· VALUES CHECKS ···

        # Ensure that the "value" is not another flag
        if value.startswith('--'):
            print(f"\nError: '{flag}' requires a value, got '{value}' instead")
            usage_instructions()

        # Ensure input_file and output_file are files and have a valid extension (.fasta or .fastq)
        if flag in ["--input", "--output"]:
            if not re.match(r'^[\w./\\ -]+\.(fasta|fastq)$', value, re.IGNORECASE):
                print(
                    f"\nError: '{value}' is not a valid file. Must be '.fasta' or '.fastq'")
                usage_instructions()

        # Ensure that the operation flag has a valid value
        if flag == "--operation":
            # Ensure that the operation is valid
            if value.lower() not in valid_operation_values:
                print(f"\nError: '{value}' is not a valid operation")
                usage_instructions()

        # Ensure that trimming lengths are integers
        if flag in ["--trim-left", "--trim-right"]:
            if not re.match(r'^\d+$', value):
                print(
                    f"\nError: '{value}' is not a valid integer for '{flag}'")
                usage_instructions()

        # Ensure that adaptor sequence is valid (only contains A, T, C, G)
        if flag == "--adaptor":
            if not re.match(r'^[ATCG]+$', value, re.IGNORECASE):
                print(
                    f"\nError: '{value}' is not a valid adaptor sequence. Only A, T, C, G are allowed.")
                usage_instructions()

        # Associate each flag with its value
        match flag:

            # input argument
            case "--input":
                input_file = value
            # output argument
            case "--output":
                output_file = value
            # operation argument
            case "--operation":
                operation = value
            # values for trimming if requested: '--operation trim'
            case "--trim-left":
                ltrim_len = int(value)
            case "--trim-right":
                rtrim_len = int(value)
            # value for adaptor sequence if requested: '--operation adaptor-removal'
            case "--adaptor":
                adaptor = value

        # Move to the next flag-value pair
        i += 2

    # Ensure input_file and output_file have the same file extension
    input_ext = input_file.lower().split('.')[-1]
    output_ext = output_file.lower().split('.')[-1]
    if input_ext != output_ext:
        print(
            f"\nError: Input file extension '.{input_ext}' does not match output file extension '.{output_ext}'")
        usage_instructions()

    return input_file, output_file, operation, ltrim_len, rtrim_len, adaptor


def check_file(in_file, out_file):
    """
    Validate input and output files and determine the input file format.

    This function performs several critical checks:
    - Verifies input file exists and is readable
    - Determines file format (FASTA or FASTQ) by checking first character
    - Checks output file permissions
    - Prompts user before overwriting existing output files

    Parameters:
        in_file (str): Path to the input FASTA or FASTQ file to be validated and processed. 
            The file must exist, be readable, and start with either '@' (FASTQ) or '>' (FASTA).
        out_file (str): Path to the output file where processed sequences will be written. 
            If the file already exists, the user will be prompted for confirmation before overwriting.

    Returns:
        str: The detected file format, either:
            - 'fastq' if the input file starts with '@'
            - 'fasta' if the input file starts with '>'
    """

    # Check input file existance and permission
    try:
        with open(in_file, 'r') as f:
            first_char = f.read(1)
    except FileNotFoundError:
        print(f"\nError: Input file '{in_file}' not found!")
        sys.exit(1)
    except PermissionError:
        print(f"\nError: Permission denied for input file '{in_file}'!")
        sys.exit(1)

    # Check output file
    try:
        with open(out_file, 'r') as f:
            print(
                f"\nWarning: Output file '{out_file}' already exists and will be overwritten!")
            while True:
                response = input(
                    "Do you want to overwrite it? (y/n): ").strip().lower()
                if response == 'y':
                    break  # Continue with overwrite
                elif response == 'n':
                    print("Operation cancelled by user.")
                    sys.exit(0)
                else:
                    print("Please enter 'y' or 'n'")
    except FileNotFoundError:
        pass  # just continues, if the output file does not exist yet (good!)

    # Check output file permission
    try:
        with open(out_file, 'w') as f:
            pass
    except PermissionError:
        print(f"\nError: Permission denied for output file '{out_file}'!")
        sys.exit(1)

    # Determine file format based on header
    if first_char == '@':
        file_format = 'fastq'
    elif first_char == '>':
        file_format = 'fasta'
    else:
        print(
            f"\nError: Input file '{in_file}' is neither FASTA nor FASTQ format!")
        sys.exit(1)

    return file_format


# ========== OPERATIONS ==========

# + + + + + + + Reverse-complement + + + + + + +

complement = {'A': 'T', 'a': 't',
              'T': 'A', 't': 'a',
              'C': 'G', 'c': 'g',
              'G': 'C', 'g': 'c',
              'N': 'N', 'n': 'n'}


# FASTQ files
def rc_fastq(in_file, out_file):
    """
    Generate reverse complement of sequences in a FASTQ file.

    This function processes a FASTQ file (4 lines per read) and generates
    the reverse complement of each sequence. Quality scores are also reversed
    to maintain correspondence with the reversed sequence.

    Parameters:
        in_file (str): Path to the input FASTQ file.
        out_file (str): Path to the output FASTQ file where 
            reverse-complemented sequences will be written.

    Returns:
        tuple: containing three elements:
            - reads_processed (int): Total number of reads (sequences) processed.
            - bases_processed (int): Total number of nucleotide bases processed across all reads.
            - base_counts (dict): Dictionary with keys 'A', 'C', 'G', 'T', 'N' containing
            the count of each base type found in the original input sequences (before complementing).
    """
    reads_processed = 0
    bases_processed = 0
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}

    with open(in_file, 'r') as infastq, open(out_file, 'w') as outfastq:
        while True:

            # Read 4 lines at a time (FASTQ format)
            header = infastq.readline().strip()  # Line 1 (@:description)
            if not header:  # End of file
                break

            sequence = infastq.readline().strip()  # Line 2 (sequence)
            plus_line = infastq.readline().strip()  # Line 3 (+)
            quality = infastq.readline().strip()  # Line 4 (quality)

            # Determine base counts
            for base in sequence.upper():
                if base in base_counts:
                    base_counts[base] += 1
                else:
                    base_counts['N'] += 1  # Handle unexpected bases

            # Reverse of sequence
            r_sequence = sequence[::-1]

            # Complement of reverse sequence
            rc_sequence = ""
            for base in r_sequence.upper():
                if base in complement:
                    rc_sequence += complement[base]
                else:
                    rc_sequence += 'N'  # Handle unexpected bases

            # Reverse of quality
            r_quality = quality[::-1]

            # Write to output FASTQ file
            outfastq.write(
                f"{header}\n{rc_sequence}\n{plus_line}\n{r_quality}\n")

            reads_processed += 1
            bases_processed += len(sequence)

    return reads_processed, bases_processed, base_counts


# FASTA files
def rc_fasta(in_file, out_file):
    """
   Generate reverse complement of sequences in a FASTA file.

    This function processes a FASTA file (2 lines per sequence) and generates
    the reverse complement of each sequence.

    Parameters:
        in_file (str): Path to the input FASTA file.
        out_file (str): Path to the output FASTA file where 
            reverse-complemented sequences will be written.

    Returns:
        tuple: containing three elements:
            - reads_processed (int): Total number of reads (sequences) processed.
            - bases_processed (int): Total number of nucleotide bases processed across all reads.
            - base_counts (dict): Dictionary with keys 'A', 'C', 'G', 'T', 'N' containing
            the count of each base type found in the original input sequences (before complementing).
    """
    reads_processed = 0
    bases_processed = 0
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}

    with open(in_file, 'r') as infasta, open(out_file, 'w') as outfasta:
        while True:

            # Read 2 lines at a time (FASTA format)
            header = infasta.readline().strip()  # Line 1 (>:description)
            if not header:  # End of file
                break
            sequence = infasta.readline().strip()  # Line 2 (sequence)

            # Determine base counts
            for base in sequence.upper():
                if base in base_counts:
                    base_counts[base] += 1
                else:
                    base_counts['N'] += 1  # Handle unexpected bases

            # Reverse of sequence
            r_sequence = sequence[::-1]

            # Complement of reverse sequence
            rc_sequence = ""
            for base in r_sequence.upper():
                if base in complement:
                    rc_sequence += complement[base]
                else:
                    rc_sequence += 'N'  # Handle unexpected bases

            # Write to output FASTA file
            outfasta.write(f"{header}\n{rc_sequence}\n")

            reads_processed += 1
            bases_processed += len(sequence)

    return reads_processed, bases_processed, base_counts


# + + + + + + + Trim + + + + + + +

# FASTQ files
def trim_fastq(in_file, out_file, left, right):
    """
    Trim specified number of bases from both ends of sequences in a FASTQ file.

    This function removes a fixed number of bases from the 5' (left) and 3' (right)
    ends of each sequence. Quality scores are trimmed accordingly to maintain
    correspondence with the trimmed sequence.

    Parameters:
        in_file (str): Path to the input FASTQ file.
        out_file (str): Path to the output FASTQ file where trimmed sequences will be written.
        left (int): Number of bases to trim from the left (5') end of each sequence.
        right (int): Number of bases to trim from the right (3') end of each sequence.

    Returns:
        tuple: containing five elements:
            - reads_processed (int): Total number of reads (sequences) processed.
            - bases_processed (int): Total number of nucleotide bases in 
            the original input sequences (before trimming).
            - base_counts (dict): Dictionary with keys 'A', 'C', 'G', 'T', 'N' containing 
            the count of each base type in the original input sequences. 
            - total_trimmed_bases (int): Total number of bases removed across all reads.
            - tbase_counts (dict): Dictionary with keys 'A', 'C', 'G', 'T', 'N' containing 
            the count of each base type in the trimmed portions.
    """
    reads_processed = 0
    bases_processed = 0
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    tbase_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    total_trimmed_bases = 0

    with open(in_file, 'r') as infastq, open(out_file, 'w') as outfastq:
        while True:

            # Read 4 lines at a time (FASTQ format)
            header = infastq.readline().strip()  # Line 1 (@:description)
            if not header:  # End of file
                break

            sequence = infastq.readline().strip()  # Line 2 (sequence)
            plus_line = infastq.readline().strip()  # Line 3 (+)
            quality = infastq.readline().strip()  # Line 4 (quality)

            # Determine base counts
            for base in sequence.upper():
                if base in base_counts:
                    base_counts[base] += 1
                else:
                    base_counts['N'] += 1  # Handle unexpected bases

            seq_length = len(sequence)

            # Ensure the total trimming length is greater than 0 and does not exceed the sequence length
            total_trim = left + right
            if total_trim == 0:
                print(f"\nError: Trimming length must be greater than 0!")
                sys.exit(1)
            elif total_trim >= seq_length:
                print(
                    f"\nError: Trimming length exceeds the length of the reads!")
                sys.exit(1)

            # Trimmed sequence & quality
            trimmed_sequence = sequence[left:seq_length - right]
            trimmed_quality = quality[left:seq_length - right]

            # Determine trimmed_portion = left_trimmed + right_trimmed
            trimmed_portions = sequence[:left] + sequence[seq_length-right:]
            # Determine trimmed base counts
            for tbase in trimmed_portions.upper():
                if tbase in tbase_counts:
                    tbase_counts[tbase] += 1
                else:
                    tbase_counts['N'] += 1

            # Write to output FASTQ file
            outfastq.write(
                f"{header}\n{trimmed_sequence}\n{plus_line}\n{trimmed_quality}\n")

            reads_processed += 1
            bases_processed += len(sequence)
            total_trimmed_bases += len(trimmed_portions)

    return reads_processed, bases_processed, base_counts, total_trimmed_bases, tbase_counts


# FASTA files
def trim_fasta(in_file, out_file, left, right):
    """
    Trim specified number of bases from both ends of sequences in a FASTA file.

    This function removes a fixed number of bases from the 5' (left) and 3' (right)
    ends of each sequence.

    Parameters:
        in_file (str): Path to the input FASTA file. 
        out_file (str): Path to the output FASTA file where trimmed sequences will be written.
        left (int): Number of bases to trim from the left (5') end of each sequence.
        right (int): Number of bases to trim from the right (3') end of each sequence.

    Returns:
        tuple: containing five elements:
            - reads_processed (int): Total number of reads (sequences) processed.
            - bases_processed (int): Total number of nucleotide bases in 
            the original input sequences (before trimming).
            - base_counts (dict): Dictionary with keys 'A', 'C', 'G', 'T', 'N' containing 
            the count of each base type in the original input sequences.
            - total_trimmed_bases (int): Total number of bases removed across all reads.
            - tbase_counts (dict): Dictionary with keys 'A', 'C', 'G', 'T', 'N' containing 
            the count of each base type in the trimmed portions.
    """
    reads_processed = 0
    bases_processed = 0
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    tbase_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    total_trimmed_bases = 0

    with open(in_file, 'r') as infasta, open(out_file, 'w') as outfasta:
        while True:

            # Read 2 lines at a time (FASTA format)
            header = infasta.readline().strip()  # Line 1 (>:description)
            if not header:  # End of file
                break
            sequence = infasta.readline().strip()  # Line 2 (sequence)

            # Determine base counts
            for base in sequence.upper():
                if base in base_counts:
                    base_counts[base] += 1
                else:
                    base_counts['N'] += 1  # Handle unexpected bases

            seq_length = len(sequence)

           # Ensure the total trimming length is greater than 0 and does not exceed the sequence length
            total_trim = left + right
            if total_trim == 0:
                print(f"\nError: Trimming length must be greater than 0!")
                sys.exit(1)
            elif total_trim >= seq_length:
                print(
                    f"\nError: Trimming length exceed the length of the reads!")
                sys.exit(1)

            trimmed_sequence = sequence[left:seq_length - right]

            # Determine trimmed_portion = left_trimmed + right_trimmed
            trimmed_portions = sequence[:left] + sequence[seq_length-right:]
            # Determine trimmed base counts
            for tbase in trimmed_portions.upper():
                if tbase in tbase_counts:
                    tbase_counts[tbase] += 1
                else:
                    tbase_counts['N'] += 1

            # Write to output FASTA file
            outfasta.write(
                f"{header}\n{trimmed_sequence}\n")

            reads_processed += 1
            bases_processed += len(sequence)
            total_trimmed_bases += len(trimmed_portions)

    return reads_processed, bases_processed, base_counts, total_trimmed_bases, tbase_counts


# + + + + + + + Adaptor-removal + + + + + + +

# FASTQ files
def adapt_removal_fastq(in_file, out_file, adaptor):
    """
    Remove adaptor sequences from the beginning of reads in a FASTQ file.

    This function searches for the specified adaptor sequence at the start of each
    read. When found, the adaptor is removed from both the sequence and its
    corresponding quality scores.

    Parameters:
        in_file (str): Path to the input FASTQ file. 
        out_file (str): Path to the output FASTQ file where processed sequences will be written.
        adaptor (str): The adaptor sequence to search for and remove from the beginning of reads. 


    Returns:
        tuple: containing four elements:
            - reads_processed (int): Total number of reads (sequences) processed.
            - bases_processed (int): Total number of nucleotide bases in 
            the original input sequences (before adaptor removal).
            - base_counts (dict): Dictionary with keys 'A', 'C', 'G', 'T', 'N' containing 
            the count of each base type in the original input sequences.
            - adaptors_removed (int): Total number of reads where the adaptor sequence 
            was found at the beginning and subsequently removed.
    """
    reads_processed = 0
    bases_processed = 0
    adaptors_removed = 0
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}

    adaptor_len = len(adaptor)
    adaptor_u = adaptor.upper()

    with open(in_file, 'r') as infastq, open(out_file, 'w') as outfastq:
        while True:

            # Read 4 lines at a time (FASTQ format)
            header = infastq.readline().strip()  # Line 1 (@:description)
            if not header:  # End of file
                break

            sequence = infastq.readline().strip()  # Line 2 (sequence)
            plus_line = infastq.readline().strip()  # Line 3 (+)
            quality = infastq.readline().strip()  # Line 4 (quality)

            # Determine base counts
            for base in sequence.upper():
                if base in base_counts:
                    base_counts[base] += 1
                else:
                    base_counts['N'] += 1  # Handle unexpected bases

            bases_processed += len(sequence)

            # Check for adaptor at the start of the sequence
            if sequence.upper().startswith(adaptor_u):
                # Remove adaptor from sequence and quality scores
                sequence = sequence[adaptor_len:]
                quality = quality[adaptor_len:]
                adaptors_removed += 1

            # Write to output FASTQ file
            outfastq.write(
                f"{header}\n{sequence}\n{plus_line}\n{quality}\n")

            reads_processed += 1

    return reads_processed, bases_processed, base_counts, adaptors_removed


# FASTA files
def adapt_removal_fasta(in_file, out_file, adaptor):
    """
   Remove adaptor sequences from the beginning of reads in a FASTA file.

    This function searches for the specified adaptor sequence at the start of each
    sequence. When found, the adaptor is removed from the sequence.

    Parameters:
        in_file (str): Path to the input FASTA file. 
        out_file (str): Path to the output FASTA file where processed sequences will be written. 
        adaptor (str): The adaptor sequence to search for and remove from the beginning of reads. 

    Returns:
        tuple containing four elements:
            - reads_processed (int): Total number of reads (sequences) processed.
            - bases_processed (int): Total number of nucleotide bases in 
            the original input sequences (before adaptor removal).
            - base_counts (dict): Dictionary with keys 'A', 'C', 'G', 'T', 'N' containing 
            the count of each base type in the original input sequences.
            - adaptors_removed (int): Total number of reads where the adaptor sequence 
            was found at the beginning and subsequently removed.
    """
    reads_processed = 0
    bases_processed = 0
    adaptors_removed = 0
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}

    adaptor_len = len(adaptor)
    adaptor_u = adaptor.upper()

    with open(in_file, 'r') as infasta, open(out_file, 'w') as outfasta:
        while True:

            # Read 2 lines at a time (FASTA format)
            header = infasta.readline().strip()  # Line 1 (>:description)
            if not header:  # End of file
                break
            sequence = infasta.readline().strip()  # Line 2 (sequence)

            # Determine base counts
            for base in sequence.upper():
                if base in base_counts:
                    base_counts[base] += 1
                else:
                    base_counts['N'] += 1  # Handle unexpected bases

            bases_processed += len(sequence)

            # Check for adaptor at the start of the sequence
            if sequence.upper().startswith(adaptor_u):
                # Remove adaptor from sequence
                sequence = sequence[adaptor_len:]
                adaptors_removed += 1

            # Write to output FASTA file
            outfasta.write(f"{header}\n{sequence}\n")

            reads_processed += 1

    return reads_processed, bases_processed, base_counts, adaptors_removed


# ============ MAIN FUNCTION ============


def main():
    """
    Main function to orchestrate FASTA/FASTQ preprocessing operations.

    This function coordinates the entire workflow:
    1. Validates command-line arguments
    2. Checks input/output files
    3. Executes the requested operation (rc, trim, or adaptor-removal)
    4. Displays comprehensive statistics about the processing

    The function handles three preprocessing operations:
        1. Reverse-complement (rc): Generate reverse complement of sequences
        2. Trim: Remove specified bases from both ends
        3. Adaptor-removal: Remove adaptor sequences from read starts

    Parameters:
        None: All parameters are obtained from command-line arguments (sys.argv).

    Returns:
        None: Results are written to the specified output file 
            and statistics are printed to standard output.
    """

    # Verify that the minimum amount of arguments are provided
    if len(sys.argv) < 7:
        usage_instructions()

    # Verify command line arguments
    input_file, output_file, operation, ltrim_len, rtrim_len, adaptor = check_arguments()

    # Verify files
    file_format = check_file(input_file, output_file)

    # ··········· --operation rc ···········
    if operation.lower() == "rc":

        # if input file is FASTQ
        if file_format == 'fastq':
            # Call the reverse-complement function for FASTQ
            reads_processed, bases_processed, base_counts = rc_fastq(
                input_file, output_file)

        # if input file is FASTA
        elif file_format == 'fasta':
            # Call the reverse-complement function for FASTA
            reads_processed, bases_processed, base_counts = rc_fasta(
                input_file, output_file)

        # Base percentages
        per_a, per_c, per_g, per_t, per_n = base_percentages(
            base_counts, bases_processed)

        # Print summary
        print(f"\n File '{input_file}' was successfully reversed-complemented ✅\n\
 Check the output file → {output_file}\n")
        print('='*65)
        print('SUMMARY'.center(65))
        print('-'*65)
        print(f" Total reads processed: {dot_number(reads_processed)}\n\
 Total bases processed: {dot_number(bases_processed)}")
        print(f" ↳({per_a}% A, {per_c}% C, {per_g}% G, {per_t}% T) | ({per_n}% N)")
        print('='*65)

    # ··········· --operation trim ···········
    if operation.lower() == "trim":

        # Ensure all necessary arguments are provided (operation trim requires at least 9 arguments)
        if len(sys.argv) < 9:
            usage_instructions()

        # if input file is FASTQ
        if file_format == 'fastq':
            # Call the trimming function for FASTQ
            reads_processed, bases_processed, base_counts, total_trimmed_bases, tbase_counts = trim_fastq(
                input_file, output_file, ltrim_len, rtrim_len)

        # if input file is FASTA
        elif file_format == 'fasta':
            # Call the trimming function for FASTA
            reads_processed, bases_processed, base_counts, total_trimmed_bases, tbase_counts = trim_fasta(
                input_file, output_file, ltrim_len, rtrim_len)

        # Base percentages
        per_a, per_c, per_g, per_t, per_n = base_percentages(
            base_counts, bases_processed)

        # Base percentages for trimmed bases
        tper_a, tper_c, tper_g, tper_t, tper_n = base_percentages(
            tbase_counts, total_trimmed_bases)

        # Print summary
        print(f"\n File '{input_file}' was successfully hard-trimmed ✅\n\
 Check the output file → {output_file}\n")
        print('='*65)
        print('SUMMARY'.center(65))
        print('-'*65)
        print(f" Total reads processed: {dot_number(reads_processed)}\n\
 Total bases processed: {dot_number(bases_processed)}")
        print(f" ↳({per_a}% A, {per_c}% C, {per_g}% G, {per_t}% T) | ({per_n}% N)")
        print('·'*65)
        print(f" Total bases trimmed: {dot_number(total_trimmed_bases)}")
        print(
            f" ↳({tper_a}% A, {tper_c}% C, {tper_g}% G, {tper_t}% T) | ({tper_n}% N)")
        print('='*65)

    # ··· ··· ··· --operation adaptor-removal ··· ··· ···
    if operation.lower() == "adaptor-removal":

        # Ensure all necessary arguments are provided (operation adaptor-removal requires 9 arguments)
        if len(sys.argv) != 9:
            usage_instructions()

        # if input file is FASTQ
        if file_format == 'fastq':
            # Call the adaptor-removal function for FASTQ
            reads_processed, bases_processed, base_counts, adaptors_removed = adapt_removal_fastq(
                input_file, output_file, adaptor)

        # if input file is FASTA
        elif file_format == 'fasta':
            # Call the adaptor-removal function for FASTA
            reads_processed, bases_processed, base_counts, adaptors_removed = adapt_removal_fasta(
                input_file, output_file, adaptor)

        # Base percentages
        per_a, per_c, per_g, per_t, per_n = base_percentages(
            base_counts, bases_processed)

        # Print summary
        print(f"\n File '{input_file}' was successfully processed ✅\n\
 Check the output file → {output_file}\n")
        print('='*65)
        print('SUMMARY'.center(65))
        print('-'*65)
        print(f" Total reads processed: {dot_number(reads_processed)}\n\
 Total bases processed: {dot_number(bases_processed)}")
        print(f" ↳({per_a}% A, {per_c}% C, {per_g}% G, {per_t}% T) | ({per_n}% N)")
        print('·'*65)
        print(f" Adaptor: {adaptor.upper()}")
        print(
            f" Total adaptors found & removed: {dot_number(adaptors_removed)}")
        print('='*65)


# Script entry point
if __name__ == "__main__":
    main()
