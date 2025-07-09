# from Bio.SeqIO.FastaIO import FastaIterator
# import matplotlib.pyplot as plt
# import os
# from Bio.SeqUtils import MeltingTemp as mt


# # File paths
# file_path = os.getcwd() + "/dengueIN.fasta"
# output_file_path = os.getcwd() + "/fastaIN.txt"

# # Lists to store GC content percentages and melting temperatures
# gc_content_list = []
# melting_temp_list = []
# fasta_data = []  # Store processed data during the first loop

# # Calculate melting temperature
# def calculate_melting_temp(gc_content,dna_length):
#     return 64.9 + 0.41 * gc_content - (500 / dna_length)
#     # return mt.Tm_GC(fasta_sequence)
# # Read fasta file and calculate GC content and melting temperature
# with open(output_file_path, 'w') as output_file:
#     with open(file_path) as handle:
#         for values in FastaIterator(handle):
#             fasta_sequence = values.seq
#             quant_a = fasta_sequence.count("A")
#             quant_c = fasta_sequence.count("C")
#             quant_t = fasta_sequence.count("T")
#             quant_g = fasta_sequence.count("G")
#             c_g = quant_g + quant_c
#             tot = quant_g + quant_c + quant_a + quant_t
#             gc_percentage = (c_g / tot) * 100
#             gc_content_list.append(gc_percentage)
#             output_file.writelines("Arquivo fasta:\n")
#             output_file.writelines(f' A: {quant_a} \n C: {quant_c} \n T: {quant_t} \n G: {quant_g} \n GC total: {c_g} \n total: {tot} percentage: {((c_g / tot) * 100):.2f}% \n')
#             fasta_data.append((fasta_sequence, gc_percentage))  # Store data for later processing
# # Process stored data and write to file
# with open(output_file_path, 'a') as output_file:
#     for i, (fasta_sequence, gc_percentage) in enumerate(fasta_data):
#         melting_temp=(calculate_melting_temp(fasta_sequence,len(fasta_sequence)))
#         melting_temp_list.append(melting_temp)
#         output_file.write(f'Seq {i}\t{gc_percentage:.2f}%\t\t{melting_temp:.2f}\n')


# # Plotting GC content vs Melting Temperature
# plt.figure(figsize=(8, 6))
# plt.scatter(gc_content_list, melting_temp_list, color='blue', alpha=0.5)
# plt.title('GC Content vs Melting Temperature')
# plt.xlabel('GC Content (%)')
# plt.ylabel('Melting Temperature (°C)')
# plt.grid(True)
# plt.show()

#algoritmo
#1- Importe os módulos necessários
#2- Defina os caminhos dos arquivos de entrada e saída
#3 - Defina uma função para calcular a temperatura de fusão com base no conteúdo GC e no comprimento do DNA
#4 - abra o arquivo fasta e itere sobre os valores e calcule a quantidade total de cada sequencia de DNA e cada base nitrogenada
#5 - calcule a porcentagem de CG
#6 - armazenar os resultados de A,C,T,G o total de nucleotideos de cada sequencia e a tabela com a temperatura de Melting de cada sequencia



from Bio.SeqIO.FastaIO import FastaIterator
import matplotlib.pyplot as plt
import os
import math

# File paths
file_path = os.getcwd() + "/dengueIN.fasta"
output_file_path = os.getcwd() + "/fastaIN.txt"

# Lists to store GC content percentages and melting temperatures
gc_content_list = []
melting_temp_list = []

# Calculate melting temperature
def calculate_melting_temp(gc_content, dna_length):
    # Constants
    Na_conc = 100  # mM
    tm = 81.5 + 16.6 * math.log10(Na_conc) + 0.41 * gc_content - (500 / dna_length)
    return tm

# Read fasta file and calculate GC content and melting temperature
with open(output_file_path, 'w') as output_file:
    with open(file_path) as handle:
        for values in FastaIterator(handle):
            fasta_sequence = values.seq
            quant_a = fasta_sequence.count("A")
            quant_c = fasta_sequence.count("C")
            quant_t = fasta_sequence.count("T")
            quant_g = fasta_sequence.count("G")
            c_g = quant_g + quant_c
            tot = quant_g + quant_c + quant_a + quant_t
            gc_percentage = (c_g / tot) * 100
            gc_content_list.append(gc_percentage)
            melting_temp_list.append(calculate_melting_temp(gc_percentage, len(fasta_sequence)))
            output_file.writelines("Arquivo fasta:\n")
            output_file.writelines(f' A: {quant_a} \n C: {quant_c} \n T: {quant_t} \n G: {quant_g} \n GC total: {c_g} \n total: {tot} percentage: {((c_g / tot) * 100):.2f}% \n')

# Plotting
plt.figure(figsize=(8, 6))
plt.scatter(gc_content_list, melting_temp_list, marker='o')
plt.title('Melting Temperature vs. GC Content')
plt.xlabel('GC Content (%)')
plt.ylabel('Melting Temperature (°C)')
plt.grid(True)
plt.show()




# # from Bio.SeqIO.FastaIO import FastaIterator
# # import matplotlib.pyplot as plt
# # import os

# # # File paths
# # file_path = os.getcwd() + "/dengueIN.fasta"
# # output_file_path = os.getcwd() + "/fastaIN.txt"

# # # Lists to store GC content percentages and melting temperatures
# # gc_content_list = []
# # melting_temp_list = []

# # # Calculate melting temperature
# # def calculate_melting_temp(gc_content, dna_length):
# #     return 64.9 + 0.41 * gc_content - (500 / dna_length)

# # # Read fasta file and calculate GC content and melting temperature
# # with open(output_file_path, 'w') as output_file:
# #     with open(file_path) as handle:
# #         for values in FastaIterator(handle):
# #             fasta_sequence = values.seq
# #             quant_a = fasta_sequence.count("A")
# #             quant_c = fasta_sequence.count("C")
# #             quant_t = fasta_sequence.count("T")
# #             quant_g = fasta_sequence.count("G")
# #             c_g = quant_g + quant_c
# #             tot = quant_g + quant_c + quant_a + quant_t
# #             gc_percentage = (c_g / tot) * 100
# #             gc_content_list.append(gc_percentage)
# #             melting_temp = calculate_melting_temp(gc_percentage, tot)
# #             melting_temp_list.append(melting_temp)
# #             output_file.writelines("Arquivo fasta:\n")
# #             output_file.writelines(f' A: {quant_a} \n C: {quant_c} \n T: {quant_t} \n G: {quant_g} \n GC total: {c_g} \n total: {tot} percentage: {gc_percentage:.2f}% \n')
# #             output_file.writelines(f'Melting Temperature: {melting_temp:.2f}°C\n')

# # # Plotting
# # plt.figure(figsize=(8, 6))
# # plt.plot(gc_content_list, melting_temp_list, marker='o', linestyle='-')
# # plt.title('Melting Temperature vs. GC Content')
# # plt.xlabel('GC Content (%)')
# # plt.ylabel('Melting Temperature (°C)')
# # plt.grid(True)
# # plt.show()
