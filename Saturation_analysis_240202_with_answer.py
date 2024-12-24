#%%
#-*- coding:utf-8 -*-
import pandas as pd
import editdistance

import os, sys

target = int(sys.argv[1])
file_name = sys.argv[2]


primer_set = pd.read_csv("/media/backup/LKS/working_space/ATM_primers_new_21_primer.txt", sep="\t", index_col=None)



FP1 = primer_set.loc[primer_set['exon'] == "{}F".format(target)].iat[0,1]
RP1 = primer_set.loc[primer_set['exon'] == "{}R".format(target)].iat[0,1]
Ref = primer_set.loc[primer_set['exon'] == "{}R".format(target)].iat[0,2]


def read_fastq(file):
    input_read = {}
    for fq_gz in file:
        with open("./{}".format(fq_gz), 'r') as f:
            for i, line in enumerate(f):
                if i%4 != 1: pass
                else:
                    x = line.strip()
                    if x in input_read.keys():
                        input_read[x] +=1
                    else:
                        input_read[x] = 1

    return input_read

def reverse_complement(sSeq):
    dict_sBases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'U': 'U', 'n': '',
                   '.': '.', '*': '*', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', '>':'>'}
    list_sSeq = list(sSeq) 
    list_sSeq = [dict_sBases[sBase] for sBase in list_sSeq]
    return ''.join(list_sSeq)[::-1]

def find_sequence(fp, rp, read):
    i = 0
    n1 = len(fp[-8:])
    n2 = len(rp[-8:])
    #let n1 < n2
    for i in range(len(read) - n1 + 1):
        subseq = read[i:i+n1]
        if editdistance.eval(fp[-8:], subseq) <= 2:
            return True
        
    rev_read = reverse_complement(read)
    for i in range(len(rev_read) - n2 + 1):
        subseq = rev_read[i:i+n2]
        if editdistance.eval(rp[-8:], subseq) <= 2:
            return True
    
    return False

def fetch_read(fp, rp, input_read_dict):
    fetched_read = {}
        
    for key, val in input_read_dict.items():
        if find_sequence(fp, rp, key) == True:
            if key in fetched_read.keys():
                fetched_read[key] = fetched_read[key] + val
            else:
                fetched_read[key] = val
        else:
            continue
    
    return fetched_read

import os
  
list_files = os.listdir('./')
files = [value for value in list_files if file_name in value and 'fastq' in value]

fastq_input = read_fastq(files)
print(sum(fastq_input.values()))
fastq_fetched = fetch_read(FP1, RP1, fastq_input)
total_reads = sum(fastq_fetched.values())

def find_best_match_location(sequence, search_string):
	
	score_list = []
	# Make list of scores for each window
	for i in range(0,len(sequence)-len(search_string)+1):
		current_string = sequence[i: i+len(search_string)]
		score = 0
		for i,base in enumerate(current_string):
			if base == search_string[i]:
				score = score
			else:
				score +=1
		score_list.append(score)

	if score_list:
		min_score_index = score_list.index(min(score_list))
	else:
		min_score_index = 0
		score_list = [None]
	
	#returns tuple of search strings starting and ending index
	#and match score
	return (min_score_index, 
			min_score_index + len(search_string), 
			min(score_list))


def remove_lowercase_letters(string):
    uppercase_string = ''.join(c for c in string if not c.islower())
    return uppercase_string
#%%
library = pd.read_csv("/media/backup/LKS/working_space/Edit_answer_with_pam.csv", index_col=None)
library_exon = library[library['Exon']=='exon_{}'.format(target)].reset_index(drop=True)

PBS_RT = {}
for idx, ser in library_exon[['Inderx','RT_PBS_sense']].iterrows():
    PBS_RT[ser['Inderx']] = ser['RT_PBS_sense']


def remove_wild_type(ref, dictionary):
    keys_to_remove = []
    a = dictionary.copy()

    for key in a.keys():
        if key in ref:
            keys_to_remove.append(key)
    for key in keys_to_remove:
        del a[key]
    
    return a

fastq_fetched_wo_WT = remove_wild_type(Ref, fastq_fetched)

#%%
def find_pbsrt(lib_seq, dictionary):
    pbsrt_cnt = []
    for pbs_rt in lib_seq['RT_PBS_sense']:
        cnt = 0
        for read, count in dictionary.items():
            if pbs_rt in read:
                cnt += count
            else: pass
        pbsrt_cnt.append(cnt)

    lib_seq['count'] = pbsrt_cnt

    return lib_seq

PE_result = find_pbsrt(library_exon, fastq_fetched_wo_WT)

#%%
def seg_generator(minimum_location, maximum_location, ref_seq):
    #ref is string sequence
    each_minimum = minimum_location-65
    each_maximum = maximum_location+65
    segment_seq = ref_seq[each_minimum : each_maximum + 1]

    return each_minimum, each_maximum, segment_seq

ATM_ref = pd.read_csv("/media/backup/LKS/working_space/ATM_gene_annotation4.txt", sep="\t", index_col=None)
ATM_ref_seq = ATM_ref.loc[ATM_ref['Name'] == 'Whole_Ref_sequence', 'sequence'].reset_index(drop = True)[0]
ATM_exon_df = ATM_ref[ATM_ref['Type'] == 'exon'].copy()
ATM_exon_df


ATM_exon_segment_df = ATM_exon_df.copy()

for idx in ATM_exon_segment_df.index:
    
    each_name = ATM_exon_df.loc[idx, 'Name']
    minimum_location = ATM_exon_df.loc[idx, 'Minimum']
    maximum_location = ATM_exon_df.loc[idx, 'Maximum']
    ref_seq = ATM_ref_seq
    
    each_minimum, each_maximum, each_exon_segment_seq = seg_generator(minimum_location, maximum_location, ref_seq)
    
    ATM_exon_segment_df.loc[idx, 'Name'] = each_name + '_seg'
    ATM_exon_segment_df.loc[idx, 'Type'] = 'exon_seg'
    ATM_exon_segment_df.loc[idx, 'Minimum'] = each_minimum
    ATM_exon_segment_df.loc[idx, 'Maximum'] = each_maximum
    ATM_exon_segment_df.loc[idx, 'sequence'] = each_exon_segment_seq
    ATM_exon_segment_df.loc[idx, 'Length'] = len(each_exon_segment_seq)

# Exon 별로 시작 부위에 코돈의 염기가 몇개 남는지 구한다 for positive strand.
left_codon_nuc_pos = [0,0]
total_nuc = 457
for i in ATM_exon_df.loc[2:,'Length']:
    left_codon_nuc_pos.append((3-(total_nuc-385)%3)%3)
    total_nuc = total_nuc + i

ATM_exon_df['left_codon_nuc_pos'] = left_codon_nuc_pos

ref_CDS = ATM_exon_df.loc[ATM_exon_df['Name']=='ATM_exon_2'].at[1, 'sequence']
ref_CDS = ref_CDS[30:]

for i in range(2,62):
    CDS = ATM_exon_df.at[i, 'sequence']
    ref_CDS = ref_CDS + CDS

CDS = ATM_exon_df.at[62, 'sequence'][:184]
ref_CDS = ref_CDS + CDS
    
codontab = {
    'TCA': 'S',    # Serine
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G',     # Glicina
    '***': '**',
    'intron': 'INTRON'
}


aminotab = {
    'TCA': ['TCC', 'TCG', 'TCT'],    # Serina
    'TCC': ['TCA', 'TCG', 'TCT'],    # Serina
    'TCG': ['TCC', 'TCA', 'TCT'],    # Serina
    'TCT': ['TCC', 'TCA', 'TCG'],    # Serina
    'TTC': ['TTT'],    # Fenilalanina
    'TTT': ['TTC'],    # Fenilalanina
    'TTA': ['TTG'],    # Leucina
    'TTG': ['TTA'],    # Leucina
    'TAC': ['TAT'],    # Tirosina
    'TAT': ['TAC'],    # Tirosina
    'TAA': ['TAG'],    # Stop
    'TAG': ['TAA'],    # Stop
    'TGC': ['TGT'],    # Cisteina
    'TGT': ['TGC'],    # Cisteina
    'TGA': ['TAA'],    # Stop
    'TGG': ['TGG'],    # Triptofano
    'CTA': ['CTC', 'CTG', 'CTT'],    # Leucina
    'CTC': ['CTA', 'CTG', 'CTT'],    # Leucina
    'CTG': ['CTA', 'CTC', 'CTT'],    # Leucina
    'CTT': ['CTA', 'CTC', 'CTG'],    # Leucina
    'CCA': ['CCC', 'CCG', 'CCT'],    # Prolina
    'CCC': ['CCA', 'CCG', 'CCT'],    # Prolina
    'CCG': ['CCA', 'CCC', 'CCT'],    # Prolina
    'CCT': ['CCA', 'CCC', 'CCG'],    # Prolina
    'CAC': ['CAT'],    # Histidina
    'CAT': ['CAC'],    # Histidina
    'CAA': ['CAG'],    # Glutamina
    'CAG': ['CAA'],    # Glutamina
    'CGA': ['CGC', 'CGG', 'CGT'],    # Arginina
    'CGC': ['CGA', 'CGG', 'CGT'],    # Arginina
    'CGG': ['CGA', 'CGC', 'CGT'],    # Arginina
    'CGT': ['CGA', 'CGC', 'CGG'],    # Arginina
    'ATA': ['ATC', 'ATT'],    # Isoleucina
    'ATC': ['ATA', 'ATT'],    # Isoleucina
    'ATT': ['ATA', 'ATC'],    # Isoleucina
    'ATG': ['ATG'],    # Methionina
    'ACA': ['ACC', 'ACG', 'ACT'],    # Treonina
    'ACC': ['ACA','ACG', 'ACT'],    # Treonina
    'ACG': ['ACA','ACC', 'ACT'],    # Treonina
    'ACT': ['ACA','ACC', 'ACG'],    # Treonina
    'AAC': ['AAT'],    # Asparagina
    'AAT': ['AAC'],    # Asparagina
    'AAA': ['AAG'],    # Lisina
    'AAG': ['AAA'],    # Lisina
    'AGC': ['AGT'],    # Serina
    'AGT': ['AGC'],    # Serina
    'AGA': ['AGG'],    # Arginina
    'AGG': ['AGA'],    # Arginina
    'GTA': ['GTC', 'GTG', 'GTT'],    # Valina
    'GTC': ['GTA', 'GTG', 'GTT'],    # Valina
    'GTG': ['GTA', 'GTC', 'GTT'],    # Valina
    'GTT': ['GTA', 'GTC', 'GTG'],    # Valina
    'GCA': ['GCC', 'GCG', 'GCT'],    # Alanina
    'GCC': ['GCA', 'GCG', 'GCT'],    # Alanina
    'GCG': ['GCA', 'GCC', 'GCT'],    # Alanina
    'GCT': ['GCA', 'GCC', 'GCG'],    # Alanina
    'GAC': ['GAT'],    # Acido Aspartico
    'GAT': ['GAC'],    # Acido Aspartico
    'GAA': ['GAG'],    # Acido Glutamico
    'GAG': ['GAA'],    # Acido Glutamico
    'GGA': ['GGC', 'GGG', 'GGT'],    # Glicina
    'GGC': ['GGA', 'GGG', 'GGT'],    # Glicina
    'GGG': ['GGA', 'GGC', 'GGT'],    # Glicina
    'GGT': ['GGA', 'GGC', 'GGG']     # Glicina
}

PE_result['vaf'] = (PE_result['count']+1)/total_reads
PE_result['rpm'] = (PE_result['count']+1)/total_reads
PE_efficiency = sum(PE_result['count'])*100/total_reads
print(PE_efficiency)
statistics = pd.DataFrame([PE_efficiency, sum(PE_result['count']), total_reads])

#%%
statistics.to_csv('./{}_stats.csv'.format(file_name), sep=',', na_rep='NaN')

PE_result.to_csv('./{}.csv'.format(file_name), sep=',', na_rep='NaN')


# Return all possible amino acids based on the codon input from the reference coding sequence

def find_possible_aa(codon):
    iterator_2 = ['A', 'T', 'G', 'C']
    changed_codon = []
    for i in range(len(codon)):
        for j in range(0,4):
            if i == 0:
                changed_codon.append(iterator_2[j] + codon[1:])
            elif i == 1:
                changed_codon.append(codon[0] + iterator_2[j] + codon[2])
            else:
                changed_codon.append(codon[:2] + iterator_2[j])

    changed_amino_acid = list(set(list(codontab[i] for i in changed_codon)))
    return changed_amino_acid

def find_possible_nuc(ref_CDS):
    iterator_2 = ['A', 'T', 'G', 'C']
    changed_nuc = []
    for i in range(len(ref_CDS)):
        for j in range(0,4):
            if iterator_2[j] == ref_CDS[i]:
                pass
            else:
                changed_nuc.append(str(i+1)+ref_CDS[i]+'>'+iterator_2[j])
    
    return changed_nuc

#%%
#possible amino acid change data frame 

original_aa = []
position_aa = []
changed_aa = []
possible_df = pd.DataFrame()

for i in range(len(ref_CDS)//3):
    changed = find_possible_aa(ref_CDS[3*i:3*i+3])
    original_aa.extend([codontab[ref_CDS[3*i:3*i+3]]] * len(changed))
    position_aa.extend([i+1]* len(changed))

    changed_aa.extend(changed)


possible_df['original'] = original_aa
possible_df['position'] = position_aa
possible_df['changed'] = changed_aa

possible_df['HGVP'] = possible_df['original'].astype(str) + possible_df['position'].astype(str) + possible_df['changed'].astype(str)

possible_df['in_library'] = possible_df['HGVP'].isin(PE_result['HGVP_intend'])

#possible_df.to_csv('./exon_{}_{}_library_aa_saturation_result.csv'.format(target,file_name), sep=',', na_rep='NaN')

# %%
possible_nuc = pd.DataFrame()
possible_nuc['HGVSc'] = find_possible_nuc(ref_CDS)
possible_nuc['in_library'] = possible_nuc['HGVSc'].isin(PE_result['HGVSc_intend'])
possible_nuc.loc[possible_nuc['in_library'] == False, 'in_library'] = possible_nuc['HGVSc'].isin(list(PE_result.loc[PE_result['HGVP_intend'] == 'other','HGVSc_intend']))
possible_nuc.loc[possible_nuc['in_library'] == False, 'in_library'] = possible_nuc['HGVSc'].isin(list(PE_result.loc[PE_result['HGVP_intend'] == 'other', 'HGVSc_syn']))

#possible_nuc.to_csv('./exon_{}_{}_library_nucleotide_saturation_result.csv'.format(target, file_name), sep=',', na_rep='NaN')