from Bio import SeqIO
from pathlib import Path
import os

#usando leitura de fasta

file_path= os.getcwd()+"/dengue.fasta"
output_file_path = os.getcwd() + "/fasta.txt"
fasta_sequence= SeqIO.read(open(file_path),"fasta")

quant_a=fasta_sequence.seq.count("A")
quant_c=fasta_sequence.seq.count("C")
quant_t=fasta_sequence.seq.count("T")
quant_g=fasta_sequence.seq.count("G")

tot= quant_g+quant_a+quant_c+quant_t

print("Arquivo fasta:")
print (f' A: {quant_a} \n C: {quant_c} \n T: {quant_t} \n G: {quant_g} \n total nucleotídeos: {tot}')

fasta_sequence= SeqIO.read(open(file_path),"fasta")

arquivo = Path("resultados.txt")

with open(output_file_path, 'w') as output_file:
    # Escrever as informações no arquivo
    output_file.write("Arquivo fasta:\n")
    output_file.write(f'A: {quant_a}\nC: {quant_c}\nT: {quant_t}\nG: {quant_g}\ntotal nucleotideos: {tot}')
    