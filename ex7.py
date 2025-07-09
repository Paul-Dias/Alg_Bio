def read_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    sequences = []
    current_seq = ""
    
    for line in lines:
        if line.startswith(">"):
            if current_seq:
                sequences.append(current_seq)
                current_seq = ""
        else:
            current_seq += line.strip()
    
    if current_seq:
        sequences.append(current_seq)
    
    return sequences

def find_overlap(s1, s2, min_length=3):
    start = 0
    
    while True:
        start = s1.find(s2[:min_length], start)
        
        if start == -1:
            return 0
        
        if s2.startswith(s1[start:]):
            return len(s1) - start
        
        start += 1

def merge_sequences(sequences):
    while len(sequences) > 1:
        max_len = -1
        best_pair = (0, 0)
        best_seq = ""
        
        for i in range(len(sequences)):
            for j in range(len(sequences)):
                if i != j:
                    overlap_len = find_overlap(sequences[i], sequences[j])
                    if overlap_len > max_len:
                        max_len = overlap_len
                        best_pair = (i, j)
                        best_seq = sequences[i] + sequences[j][overlap_len:]
        
        i, j = best_pair
        sequences[i] = best_seq
        del sequences[j]
    
    return sequences[0]

def write_fasta(sequence, file_path):
    with open(file_path, 'w') as file:
        file.write(">Contig\n")
        file.write(sequence + "\n")

# Path to the input file
reads_file_path = 'reads.fasta'

# Read the sequences from the input file
sequences = read_fasta(reads_file_path)

# Merge the sequences to form the contig
contig = merge_sequences(sequences)

# Path to the output file
contig_file_path = 'contig.fasta'

# Write the contig to the output file in FASTA format
write_fasta(contig, contig_file_path)
