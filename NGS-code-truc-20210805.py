#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""

Created on Thu Aug  5 15:02:30 2021



@author: trucdogunther

"""



import os



from Bio import SeqIO



file_dir = "/Users/alina/.spyder-py3"



os.chdir(file_dir)





WT_sequence_aa = "GDGSIAEYFNRPMHIHDWQIMDSGYYG"



# Do NOT include stop codon in the DNA sequence

WT_sequence_DNA = "GGCGATGGCAGCATTGCGGAATACTTTAACCGTCCGATGCATATTCATGATTGGCAGATTATGGATAGCGGCTATTATGGC"





file_name = "AT-11-DecNGS"



# open and parse the input FASTA file 

fasta_file = SeqIO.parse(file_name+".fasta", "fasta")



sequence_list = []



for sequence in fasta_file:

    

    # number of times current sequence was sequenced

    count = int(sequence.id.split("-")[1])

    # amino acid sequence of peptide

    peptide_sequence = str(sequence.translate().seq)

    DNA_sequence = str(sequence.seq)

    

    sequence_list.append((count , peptide_sequence , DNA_sequence))

    



# key is WT amino acid-amino acid position-mutant amino acid

# value is int of number of time mutation(s) appears

mutations_dict_aa = {}



# key is WT codon-start nucleotide position-mutant condon

# value is tuple of (int of number of time codon mutant appears , WT codon, str of nucleotide position, mutant codon)

mutations_dict_DNA = {}

    

for sequence_tuple in sequence_list:

    

    # compare amino acid sequence first

    

    if len(sequence_tuple[1]) == len(WT_sequence_aa):

    

        aa_counter = 0

        

        WT_flag_aa = True

        

        mutation_tracker_aa = ""

        

        for amino_acid in WT_sequence_aa:

            

            if amino_acid != sequence_tuple[1][aa_counter]:

                

                mutation_tracker_aa = mutation_tracker_aa + amino_acid + str(aa_counter+2) + sequence_tuple[1][aa_counter]

            

                WT_flag_aa = False

                

            aa_counter = aa_counter + 1

            

        mutation_count_aa = mutations_dict_aa.get(mutation_tracker_aa,0)

                

        if WT_flag_aa == False:

            mutations_dict_aa[mutation_tracker_aa] = mutation_count_aa + sequence_tuple[0]

        else:  

            WT_count_aa = mutations_dict_aa.get("WT",0)

            

            mutations_dict_aa["WT"] = WT_count_aa + sequence_tuple[0]

        

    

    # compare DNA sequence next

            

    if len(sequence_tuple[2]) == len(WT_sequence_DNA):

                    

        for nucleotide_position in range(0, len(WT_sequence_DNA), 3):

            

            WT_codon = WT_sequence_DNA[nucleotide_position : nucleotide_position+3]

            

            library_codon = sequence_tuple[2][nucleotide_position : nucleotide_position+3]

            

            if WT_codon != library_codon:

                

                mutation_DNA = WT_codon + str(nucleotide_position + 1) + library_codon

                

                if mutations_dict_DNA.get(mutation_DNA , 0) == 0:

                    mutation_count_DNA = 0

                else:

                    mutation_count_DNA = mutations_dict_DNA[mutation_DNA][0]

                                                

                mutations_dict_DNA[mutation_DNA] = (mutation_count_DNA + sequence_tuple[0] , WT_codon, str(nucleotide_position + 1), library_codon)

                

 

# print out amino acid mutation file

output_file_name_1 = file_name + "_output_aminoacid.txt"

os.chdir(file_dir)

output_file_handle_1 = open(output_file_name_1, "w")



output_file_handle_1.write("Mutation\tMutation count\n")



for mutation_aa,count_aa in mutations_dict_aa.items():

    

    output_file_handle_1.write(mutation_aa + "\t")

    output_file_handle_1.write(str(count_aa) + "\n")

    

output_file_handle_1.close()



# print out DNA sequence mutation file

output_file_name_2 = file_name + "_output_DNA.txt"

os.chdir(file_dir)

output_file_handle_2 = open(output_file_name_2, "w")



output_file_handle_2.write("Codon mutation\tWT codon\t Start nucleotide position\tMutant codon\tMutation count\n")



for mutation_DNA,mutant_tuple in mutations_dict_DNA.items():

    

    output_file_handle_2.write(mutation_DNA + "\t")

    output_file_handle_2.write(mutant_tuple[1] + "\t")

    output_file_handle_2.write(mutant_tuple[2] + "\t")

    output_file_handle_2.write(mutant_tuple[3] + "\t")

    output_file_handle_2.write(str(mutant_tuple[0]) + "\n")

    

output_file_handle_2.close()

    