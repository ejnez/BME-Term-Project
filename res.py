import argparse
import math
import pandas as pd
from Bio import SeqIO
import re

vals = 'ACDEFGHIKLMNPQRSTVWY-'

def fasta_to_array(fasta_file):
    """Converts the sequences in an aligned FASTA file to a list of (formatted_species_name, sequence) tuples
    using the full species name from the 'OS' field in the FASTA header and formats it.
    """
    sequences = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract the species name from the OS field in the description
        description = record.description
        os_start = description.find("OS=")
        if os_start != -1:
            os_end = description.find("OX=", os_start)  # Find end of OS field (before OX field)
            full_species_name = description[os_start + 3:os_end].strip()
            
            # Split the species name into parts
            name_parts = full_species_name.split()
            
            # Format the species name as G.species, e.g., C.sabaeus
            if len(name_parts) >= 2:
                genus_initial = name_parts[0][0]  # First letter of the genus
                species_part = ''.join(name_parts[1:])  # Combine species and subspecies names if present
                formatted_species_name = f"{genus_initial}.{species_part}"
            else:
                formatted_species_name = full_species_name  # Fallback to full name if splitting fails
        else:
            formatted_species_name = record.id  # Fallback to using sequence ID if OS field is not found

        # Store the formatted species name and its corresponding sequence
        sequences.append((formatted_species_name, str(record.seq)))

    return sequences

def calculate_pwm(sequences, pseudocount=1):
    """Calculates a position weight matrix from array of aligned 
    sequences.
    """
    counts = {residue: [0] * len(sequences[0]) for residue in vals}

    for sequence in sequences:
        for i, residue in enumerate(sequence):
            if residue in counts:
                counts[residue][i] += 1
            else:
                print(f"Invalid character '{residue}' in sequence. Skipping.")

    for residue in counts:
        for i in range(len(counts[residue])):
            counts[residue][i] += pseudocount

    pwm = {residue: [count / (len(sequences) + pseudocount * len(
        counts)) for count in counts[residue]] for residue in counts}
    return pwm

def calculate_background_frequency(sequences):
    """Calculates background frequency of each residue (or gap)
    over the aligned sequences.
    """
    dictbg = {val: 0 for val in vals}
    dictbg["total"] = 0

    for line in sequences:
        for i in line:
            if i in dictbg:
                dictbg[i] += 1
                dictbg["total"] += 1
            else:
                print(f"Invalid character '{i}' in sequence. Skipping.")

    dictf = {i: dictbg[i] / dictbg["total"] for i in vals}
    return dictf

def calculate_position_score(
        sequence, pwm, background_freq, gap_open_penalty=2.90, gap_extend_penalty=0):
    """Calculates the position-specific score for a sequence 
    using pwm and background frequency.
    """
    num_positions = len(sequence)
    scores = []

    gap_opened = False

    for i in range(num_positions):
        score = 0.0

        residue = sequence[i]
        if residue == '-':
            if not gap_opened:
                score -= gap_open_penalty
                gap_opened = True
            else:
                score -= gap_extend_penalty
        else:
            gap_opened = False
            if residue in pwm and i < len(pwm[residue]) and residue in background_freq:
                if pwm[residue][i] > 0 and background_freq[residue] > 0:
                    score += math.log2(pwm[residue][i] / background_freq[residue])

        scores.append(score)
        print(f"Score for position {i + 1}: {score}")

    return scores

def extract_os_name(header):
    """Extract the OS (organism name) from the FASTA header."""
    match = re.search(r'OS=([^ ]+)', header)  # Extract OS field
    if match:
        return match.group(1)
    else:
        return header  # If OS not found, return the full header

def calculate_and_save_scores_to_csv(sequences_with_headers, output_file):
    """Calculate and write scores to a CSV file, excluding positions with gaps in more than 5 sequences.
    
    Each column represents scores for a sequence, and rows contain residue-wise scores 
    with labels for each position in the aligned sequence.
    """

    headers = [header for header, _ in sequences_with_headers]
    sequences = [sequence for _, sequence in sequences_with_headers]

    pwm = calculate_pwm(sequences)
    background_freq = calculate_background_frequency(sequences)

    all_scores = []

    for sequence in sequences:
        sequence_scores = calculate_position_score(sequence, pwm, background_freq)
        all_scores.append(sequence_scores)

    df = pd.DataFrame(all_scores).transpose()

    os_names = [extract_os_name(header) for header in headers]
    df.columns = [f'{os_name}' for os_name in os_names]

    seq_df = pd.DataFrame([list(seq) for seq in sequences])

    gaps_per_position = (seq_df == '-').sum(axis=0)
    positions_to_keep = gaps_per_position[gaps_per_position <= 5].index  

    df_filtered = df.loc[positions_to_keep]
    varying_positions = df_filtered.apply(lambda row: row.nunique() > 1, axis=1)  
    varying_scores_df = df_filtered[varying_positions].copy()

    varying_scores_df['Position'] = df_filtered.index[varying_positions] + 1

    cols = ['Position'] + [col for col in varying_scores_df.columns if col != 'Position']
    varying_scores_df = varying_scores_df[cols]

    transposed_df = varying_scores_df.set_index('Position').transpose()

    transposed_df.to_csv(output_file)

    print(f"Scoring completed and saved to {output_file}.")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate sequence scores and output to CSV.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file with aligned sequences.')
    parser.add_argument('output_file', type=str, help='Output CSV file for scores.')
    args = parser.parse_args()

    sequences_with_headers = fasta_to_array(args.fasta_file)
    
    calculate_and_save_scores_to_csv(sequences_with_headers, args.output_file)

