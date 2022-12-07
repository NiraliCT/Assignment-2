import os, sys, random
from urllib.request import urlopen


#data extract from URL to dictionary
fasta_dict = {}
url = "https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/MonkeyPox.fn"
file = urlopen(url)
with urlopen(url) as f:
    for line in file:
        decoded_line = line.decode("utf-8")
        decoded_line = decoded_line.strip()
        if decoded_line == "":
            continue
        if decoded_line[0] == ">":
            fasta_dict[decoded_line] = ""
            current = decoded_line
        else:
            if fasta_dict[current] == "":
                fasta_dict[current] = decoded_line
            else:
                fasta_dict[current] += decoded_line
'''print(len(fasta_dict.keys()))'''


#Create subsequences dictionary and add random mutation and reverse complement the dna sequences
subseq_fasta = {}
for i,(k,v) in enumerate(fasta_dict.items()):
    random_key, random_value = random.choice(list(fasta_dict.items()))
    #print(len(random_value))
    random_mutation = random.choice(list('ATGC'))
    pos = random.randint(0, len(random_value) - 1)
    mut_random_value = "".join((random_value[:pos], random_mutation, random_value[pos:]))
    #print(len(mut_random_value))
    subseq_fasta[random_key] = (mut_random_value[::-1].replace("A","t").replace("T","a").replace("G","c").replace("C","g").upper())
    random_key, random_value, mut_random_value = '' , '' , ''
    i+=1
    if i == 3:
        break
with open('subsequence_fasta.fasta','w') as file:
    for key_1,Value_1 in subseq_fasta.items():
        file.write(str(key_1 + '\n' + Value_1 + '\n'))


#Create bed file
with open('subsequence_bed.bed','w') as file1:
    for indexer,name in enumerate(subseq_fasta.keys()):
        file1.write(str(name + ' : '))
        file1.write(str(indexer))
        file1.write(str(' : '))
        new_index = indexer + 100
        file1.write(str(new_index))
        file1.write(str('\n'))


#Remove subsequences from the original file
for k,v in list(fasta_dict.items()):
        if k in subseq_fasta.keys():
            del fasta_dict[k]
'''print(len(fasta_dict.keys()))'''


#Add mutated data + data from which substring data was removed + write into final out.fasta file
merged_data = {**fasta_dict, **subseq_fasta}
with open('out.fasta','w') as fhl1:
    for k,v in merged_data.items():
        fhl1.write(str(k + '\n' + v + '\n'))
#print(len(merged_data.keys()))
print('Assignment 2 Completed')
