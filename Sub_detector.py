import sys
import os
import regex
import pandas as pd
import openpyxl as op
import numpy as np

from collections import Counter
import multiprocessing

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio.SeqIO import write
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import reverse_complement

from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl import load_workbook
from openpyxl.styles import Font, Alignment





ind = sys.argv[1]
os.chdir('/home/baelab/HW_py/1_4_Subdetector/' + '/' + ind)
f_list = open('target_list.txt').readlines()

first_line_list = f_list[0][:-1].split('\t')

minumum_num=int(first_line_list[1])
comparison_range=int(first_line_list[3])
indicator_range=int(first_line_list[5])
window_range=int(first_line_list[7])







# minimum_frequency : Time Complexity O(1)

def minimum_frequency(records, n=2):
    seq_form_fastq = [record[1] for record in records]  
    counts = Counter(seq_form_fastq)    
    minimum_records = [record for record in records if counts[record[1]] >= n] 
    return minimum_records 



# indicator selection : w/ mismatch

def pick_indicator_with_mismatch(seq_input, left_indicator, right_indicator): # seq_input type = minimum_records
    indicator_records=[]
    if spacer_seq in wt_ref_seq:
        for seq in seq_input:
            is_WT_L = regex.search(rf"({left_indicator}){{s<={2}}}", seq[1]) # regex fuzzy matching LI
            is_WT_R = regex.search(rf"({right_indicator}){{s<={2}}}", seq[1]) # regex fuzzy matching RI
            if is_WT_L is not None and is_WT_R is not None:
                nomean=['A','B','C']
                nomean[1]=seq[1][is_WT_L.start():is_WT_R.end()]
                indicator_records.append(nomean)
                
    elif spacer_seq in wt_ref_seq_complement:
        for seq in seq_input:
            is_WT_L = regex.search(rf"({left_indicator}){{s<={2}}}", str(Seq(seq[1]).reverse_complement())) # regex fuzzy matching LI
            is_WT_R = regex.search(rf"({right_indicator}){{s<={2}}}", str(Seq(seq[1]).reverse_complement())) # regex fuzzy matching RI
            if is_WT_L is not None and is_WT_R is not None:
                nomean=['A','B','C']
                nomean[1]=str(Seq(seq[1]).reverse_complement())[is_WT_L.start():is_WT_R.end()]
                indicator_records.append(nomean)

    return indicator_records 






# count mutagen

def count_mutagen(seq_input, ref_seq): # seq_input type = indicator_records
    sigle_read_muts = {f'{i} mut': 0 for i in range(1, 11)}
    sub_rate = [[0,0,0,0, ref_seq[i]] for i in range(len(ref_seq))] #ATGC
    cnt = {'dimer' : 0, 'ins' : 0, 'del' : 0, 'WT' :0, 'mutants': 0, 'mutation(nucleotide)_number' : 0}
    
    
    for seq in seq_input: 
        if len(seq[1]) < len(ref_seq)*0.7:
            cnt['dimer'] += 1

        else:
            alignments = pairwise2.align.globalms(ref_seq, seq[1], 1, -1, -10, -.5)
            split_format = format_alignment(*alignments[0]).split('\n')
            alignments_list = [s.strip() for s in split_format]

            if '-' in alignments_list[0]:
                cnt['ins'] += 1
                
            elif '-' in alignments_list[2]:
                cnt['del'] += 1
                
            elif '.' not in alignments_list[1]:
                cnt['WT'] += 1
                
            else:
                cnt['mutants'] += 1
                num_mutations = alignments_list[1].count(".")
                cnt['mutation(nucleotide)_number'] += num_mutations
                mut_positions = tuple((index) for index, char in enumerate(alignments_list[1]) if char == '.')
                if num_mutations >= 10:
                    sigle_read_muts['10 mut'] += 1
                else:
                    sigle_read_muts[f'{num_mutations} mut'] += 1
                
                for i in mut_positions:
                    if alignments_list[2][i] == 'A':
                        sub_rate[i][0] += 1
                    elif alignments_list[2][i] == 'T':
                        sub_rate[i][1] += 1
                    elif alignments_list[2][i] == 'G':
                        sub_rate[i][2] += 1
                    elif alignments_list[2][i] == 'C':
                        sub_rate[i][3] += 1
    return sub_rate, cnt, sigle_read_muts
    
   








### 기록 시작부분 ###








wb = Workbook()

thisis2=2
thisis1=1

with pd.ExcelWriter('sub.xlsx', engine='openpyxl', mode='w') as writer:
    writer.book = wb
    
    for t in f_list[1:]:    
        each_line_list = t.split('\t')

        if len(each_line_list) <= 1:
            continue









################################## 기본 변수 설정 ##################################

        

        wt_ref = each_line_list[0].strip()
        spacer = each_line_list[1].strip()
        direc = each_line_list[2].strip()
        inedx = each_line_list[3].strip()
        name = each_line_list[4].strip()


        wt_ref_seq = wt_ref.upper()
        spacer_seq = spacer.upper()
        wt_ref_seq_complement=str(Seq(wt_ref_seq).reverse_complement())
        spacer_seq_complement=str(Seq(spacer_seq).reverse_complement())



        if spacer_seq in wt_ref_seq:
            position_nick = wt_ref_seq.find(spacer_seq) + 17
            
            
            if position_nick-comparison_range >= 0:
                left_indicator = wt_ref_seq[position_nick-comparison_range : position_nick-comparison_range+indicator_range]
            else:
                left_indicator = wt_ref_seq[:indicator_range]

            if position_nick+comparison_range < len(wt_ref_seq) :
                right_indicator = wt_ref_seq[position_nick+comparison_range-indicator_range : position_nick+comparison_range]
            else:
                right_indicator = wt_ref_seq[-indicator_range:]
            
            left_index = wt_ref_seq.find(left_indicator)
            right_index = wt_ref_seq.find(right_indicator)

            wt_ref_seq = wt_ref_seq[left_index:right_index+len(right_indicator)]
                
            position_nick = wt_ref_seq.find(spacer_seq) + 17
            A_list = []
            A_index=position_nick+3
            
            if position_nick+3-window_range >= 0:
                for i in range(position_nick+3-window_range, position_nick+3):  
                    if wt_ref_seq[i] == 'A':
                        A_list.append(i + 1)
            else:
                for i in range(0, position_nick+3):  
                    if wt_ref_seq[i] == 'A':
                        A_list.append(i + 1)         
            
            
            
            
            
            
        elif spacer_seq in wt_ref_seq_complement:
            position_nick = wt_ref_seq_complement.find(spacer_seq) + 17

            if position_nick-comparison_range >= 0:
                left_indicator = wt_ref_seq_complement[position_nick-comparison_range : position_nick-comparison_range+indicator_range]
            else:
                left_indicator = wt_ref_seq_complement[:indicator_range]

            if position_nick+comparison_range < len(wt_ref_seq) :
                right_indicator = wt_ref_seq_complement[position_nick+comparison_range-indicator_range : position_nick+comparison_range]
            else:
                right_indicator = wt_ref_seq_complement[-indicator_range:]
                
                
            left_index = wt_ref_seq_complement.find(left_indicator)
            right_index = wt_ref_seq_complement.find(right_indicator)

            wt_ref_seq_complement= wt_ref_seq_complement[left_index:right_index+len(right_indicator)]

            position_nick = wt_ref_seq_complement.find(spacer_seq) + 17
        
            
            A_list = []
            A_index=position_nick+3
            
            if position_nick+3-window_range >= 0:
                for i in range(position_nick+3-window_range, position_nick+3):  
                    if wt_ref_seq_complement[i] == 'A':
                        A_list.append(i + 1)
            else:
                for i in range(0, position_nick+3):  
                    if wt_ref_seq_complement[i] == 'A':
                        A_list.append(i + 1)



        else:
            print("No spacer found.")
            print(spacer_seq)
            print(wt_ref_seq_complement)
            print(f"\u25B6 \u25B6 # {inedx} analysis failed.")
            print("--------------------------------------------------------------------------")

            continue





################################## indicator를 기반으로 다시 Ref 슬라이싱 및 nick_pos 재정립 ##################################


            



################################## fastqjoin 파일 탐색 후 open ##################################


        file_name = str(inedx) + '.fastqjoin'
        os.chdir('/home/baelab/MGEL/' + direc + '/fastqjoin')


################################## fastqjoin 파일 탐색 후 open ##################################










        records = [] 
        with open(file_name, "r") as handle:
            for record in FastqGeneralIterator(handle): 
                records.append(record)



        os.chdir('/home/baelab/HW_py/1_4_Subdetector/' + '/' + ind)


        record_minimum=minimum_frequency(records, n=minumum_num)
        record_indicator=pick_indicator_with_mismatch(record_minimum, left_indicator, right_indicator)
        print(f"#{inedx} : indicator matched = {len(record_indicator)}")

        if len(record_indicator) == 0:
            print(f"\u25B6 # {inedx} has no read.")
            print("--------------------------------------------------------------------------")
            continue
        
        
        
        
        
        

        else:
            if spacer_seq in wt_ref_seq:
                record_mutagen=count_mutagen(record_indicator, wt_ref_seq)

                
            
            else:
                record_mutagen=count_mutagen(record_indicator, wt_ref_seq_complement)

                            
            
                         


            



        
        
                          
            sub_rate=record_mutagen[0]
            cnt=record_mutagen[1]
            sigle_read_muts=record_mutagen[2]
            
            mutation_number = cnt['mutants']+cnt['WT']
            indel_total = cnt['mutants']+cnt['WT']+cnt['ins']+cnt['del']
            indelper = ((cnt['ins']+cnt['del']) / indel_total * 100)
            
            if mutation_number != 0:
                A_sub_rate_percentage = [[ sub_rate[i][0] / mutation_number * 100 ] for i in range(len(sub_rate))]
                T_sub_rate_percentage = [[ sub_rate[i][1] / mutation_number * 100 ] for i in range(len(sub_rate))]
                G_sub_rate_percentage = [[ sub_rate[i][2] / mutation_number * 100 ] for i in range(len(sub_rate))]
                C_sub_rate_percentage = [[ sub_rate[i][3] / mutation_number * 100 ] for i in range(len(sub_rate))]
                sub_rate_percentage = [[ sum(sub_rate[i][:3]) / mutation_number * 100 ] for i in range(len(sub_rate))]

            else:
                print(mutation_number)
                continue

                
            ##data0 = {
            ##    'title': list(cnt.keys()),
            ##    'contents': list(cnt.values())
            ##}




            if spacer_seq in wt_ref_seq:

                data = {
                    'Seq': [item[0] for item in wt_ref_seq],
                    'total(sum)': [item[0] for item in sub_rate_percentage],
                    'A rates': [item[0] for item in A_sub_rate_percentage],
                    'T rates': [item[0] for item in T_sub_rate_percentage],
                    'G rates': [item[0] for item in G_sub_rate_percentage],
                    'C rates': [item[0] for item in C_sub_rate_percentage],

                }                
            
            else:

                data = {
                    'Seq': [item[0] for item in wt_ref_seq_complement],
                    'total(sum)': [item[0] for item in sub_rate_percentage],
                    'A rates': [item[0] for item in A_sub_rate_percentage],
                    'T rates': [item[0] for item in T_sub_rate_percentage],
                    'G rates': [item[0] for item in G_sub_rate_percentage],
                    'C rates': [item[0] for item in C_sub_rate_percentage]

                }



            ## original count mutagen (A window only)

            ##data = {
            ##    'Seq': [f"A{abs(i - A_index-1)}" for i in A_list],
            ##    'A rates': [item[0] for item in A_sub_rate_percentage],
            ##    'T rates': [item[0] for item in T_sub_rate_percentage],
            ##    'G rates': [item[0] for item in G_sub_rate_percentage],
            ##    'C rates': [item[0] for item in C_sub_rate_percentage],
            ##    'total(sum)': [item[0] for item in sub_rate_percentage],
            ##}





            ##data1 = {
            ##    'title': list(sigle_read_muts.keys()),
            ##    'contents': list(sigle_read_muts.values())
            ##}


            header0 = {
                'Min.Freq_name': ['Min.Freq'],
                'wo_name' : ['w/o dimer'],
                'indel_name': ['indel%'],
                'Memo_name': ['Memo'],
                'index_name' : ['index']
            }




            header = {
                'Min.Freq.': [str(len(record_indicator))],
                'wo': [str(mutation_number)],
                'indel': [str(round(indelper,2))],
                'file name': [str(name)],
                'index' : [str(inedx)]
            }


            #df0 = pd.DataFrame(data0)
            df = pd.DataFrame(data)
            df_header=pd.DataFrame(header)
            df_header0=pd.DataFrame(header0)
            #df1 = pd.DataFrame(data1) 

            df['A rates'] = df['A rates'].round(2)
            df['T rates'] = df['T rates'].round(2)
            df['G rates'] = df['G rates'].round(2)
            df['C rates'] = df['C rates'].round(2)
            df['total(sum)'] = df['total(sum)'].round(2)


            # df0_transposed = df0.transpose()
            df_transposed = df.transpose()
            # df1_transposed = df1.transpose()





            # base_index = {'A': 0, 'T': 1, 'G': 2, 'C': 3}

            # # 4x4 matrix
            # matrix = np.zeros((4, 4), dtype=float)

            # # sub_rate
            # for i in range(len(sub_rate)):
            #     base_from = sub_rate[i][4]
            #     for j in range(4):
            #         base_to = 'ATGC'[j]
            #         matrix[base_index[base_from]][j] += sub_rate[i][j]

            # # mutation_number
            # if mutation_number != 0:
            #     matrix_percentage = matrix / cnt['mutation(nucleotide)_number'] * 100
            # else:
            #     matrix_percentage = matrix

            # # Pandas DataFrame
            # columns = ['A', 'T', 'G', 'C']
            # index = ['\u25B6A', '\u25B6T', '\u25B6G', '\u25B6C']

            # df_matrix = pd.DataFrame(matrix_percentage.T, columns=columns, index=index).round(2)
            # df_matrix_transposed=df_matrix.transpose()





            # ?  ? ? 일??? 기
            
            # df0_transposed.to_excel(writer, sheet_name=file_name, startrow=8, startcol=0, header=False)
            df_transposed.to_excel(writer, startrow=thisis1, startcol=5, header=False)
            df_header0.to_excel(writer, startrow=1, startcol=0, header=False, index=False)
            df_header.to_excel(writer, startrow=thisis2, startcol=0, header=False, index=False)
            # df1_transposed.to_excel(writer, sheet_name=file_name, startrow=11, startcol=0, header=False)
            # df_matrix_transposed.to_excel(writer, sheet_name=file_name, startrow=2, startcol=0, header=True)
            

            thisis2 += 6
            thisis1 += 6
            
            print(f"\u25B6 #{file_name} analysis completed.")
            print("--------------------------------------------------------------------------")

print('Jobs done! mady by HW')

