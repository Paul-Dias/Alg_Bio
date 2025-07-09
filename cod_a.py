from Bio import SeqIO
import os

output_file_path=os.getcwd() + "/dengue.fasta"
output_file_path_2=os.getcwd() + "/txt.txt"

with open(output_file_path, 'r') as file:
    # Escrever as informações no arquivo
   data = file.readlines()
   print(data)
   print(data.count("A"))

count_a=0
count_c=0
count_t=0
count_g=0
for index,linha in enumerate(data):
   if index != 0:
      count_a+= linha.count("A")
      count_c+= linha.count("C")
      count_t+= linha.count("T")
      count_g+= linha.count("G")

tot=count_a+count_c+count_t+count_g

with open(output_file_path_2, 'w') as output_file:
    # Escrever as informações no arquivo
    output_file.write("Arquivo normal:\n")
    output_file.write(f'A: {count_a}\n C: {count_c}\n T: {count_t}\n G: {count_g}\ntotal nucleotideos: {tot}')