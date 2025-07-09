## Algoritmo :

# Criar funções de tradução: defina duas funções principais para traduzir sequências de RNA em sequências de proteínas. A função translate_rna_to_protein traduz uma única sequência de RNA em uma sequência de proteína, enquanto translate_rna_to_protein_all_frames traduz todas as três molduras de leitura possíveis de uma sequência de RNA em sequências de proteínas correspondentes.

# Listar e armazenar sequências de proteínas: criar uma lista para armazenar as sequências de proteínas traduzidas em diferentes molduras de leitura.

# Abertura dos arquivos de entrada e saída: abrir o arquivo de saída no modo de escrita e o arquivo de entrada FASTA contendo as sequências de DNA.

# Iteração sobre as sequências de DNA: iterar sobre as sequências de DNA no arquivo de entrada usando um iterador FASTA. Para cada sequência de DNA, é realizada a transcrição para RNA.

# Tradução e gravação das sequências de proteínas: traduzir o RNA transcrita em sequências de proteínas em todas as três molduras de leitura. Em seguida, para cada sequência de proteína em cada moldura de leitura, identifica todas as subsequências de proteínas começando com um "M" (códon de início) e terminando com um "_" (códon de parada). Cada subsequência de proteína é gravada no arquivo de saída no formato FASTA, juntamente com informações sobre a sequência original, moldura de leitura e número da proteína.

# Fechamento dos arquivos: Finalmente, os arquivos de entrada e saída são fechados.


from Bio.SeqIO.FastaIO import FastaIterator
import os
from Bio.Seq import Seq

# File paths
file_path = os.getcwd() + "/dengue_ex5.fasta"
output_file_path = os.getcwd() + "/sequencias.fasta"

# Dicionário de codificação de aminoácidos baseado em RNA
rna_codon_table = {
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
    'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',
}

# Função para traduzir RNA em proteína considerando todos os frames
def translate_rna_to_protein_all_frames(rna_sequence):
    protein_sequences = []
    
    for frame in range(3):
        protein_sequence = translate_rna_to_protein(rna_sequence[frame:])
        protein_sequences.append(protein_sequence)
        
    return protein_sequences


# Função para traduzir RNA em proteína
def translate_rna_to_protein(rna_sequence):
    protein_sequence = ""
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        if len(codon) == 3:
            protein_sequence += rna_codon_table.get(codon.upper(), 'X')
    return protein_sequence


# Listas para armazenar sequências de aminoácidos em diferentes frames
protein_sequences_frames = [[] for _ in range(3)] 

with open(output_file_path, 'w', encoding='utf-8') as output_file:
    with open(file_path) as handle:
        for values in FastaIterator(handle):
            fasta_dna = Seq(values.seq).transcribe()     
            protein_sequences = translate_rna_to_protein_all_frames(fasta_dna)
            for frame, seq in enumerate(protein_sequences):
                start = 0
                idx = 0
                while start < len(seq):
                    start = seq.find('M', start)
                    if start == -1:
                        break  # Não foi encontrada mais "M"
                    end = seq.find('_', start)
                    if end == -1:
                        end = len(seq)  # Se não houver "_", consideramos o fim da sequência
                    if start != -1:
                        segment = seq[start:end+1]
                        output_file.write(f"> {values.id}, Frame {frame+1}, Proteína {idx+1}\n")
                        output_file.write(f"{segment}\n")
                        idx += 1
                    start = end + 1  # Avançar o início para buscar a próxima sequência

