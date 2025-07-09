from Bio import SeqIO

fasta_file_path ="C:\\Users\\ppaul\\OneDrive\\Área de Trabalho\\python\\Alg_Bio\\dengue_ex5.fasta"
fasta_sequence = SeqIO.read(fasta_file_path, "fasta")

sequence = fasta_sequence.seq

translation = sequence.translate()
print(translation)

protein_lists = []
for frame in range(3):
    length = len(translation[frame:])

    for i in range(length):
        if translation[frame:][i] == "M":
            protein = ""

            for j in range(i, length):
                protein += translation[frame:][j]

                if translation[frame:][j + 1] == "*":
                    protein_lists.append(protein)
                    break

print(protein_lists)
