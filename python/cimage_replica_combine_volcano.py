import argparse
from collections import defaultdict
from scipy import stats
import scipy.stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import os
import math

###### arguments ######
parser = argparse.ArgumentParser()
parser.add_argument('-b','--basename', metavar='Str',dest='basename',help='basename',type=str,required=True)
parser.add_argument('-q','--query', metavar='Str', dest='query',help='heavy or medium',type=str, required=True)
args = parser.parse_args()
###### arguments ######

def plot_volcano(basename, query):

    starttime = datetime.datetime.now()
    print('Calculating now...')
    print('This script only supports SITES results now...')
    print('This script only supports SITES results now...')
    print('This script only supports SITES results now...')

    df = combination_filter(basename, query)
    f, ax = plt.subplots(figsize=(7, 7))

    group_num = len(df.columns) - 5
    i = 1
    a = []
    while i < group_num + 1:
        a.append('ratio_{}'.format(i))
        i += 1

    x_axis_value = df.loc[:, 2+group_num]
    y_axis_value = df.loc[:, 4+group_num]
    noinf = []
    for i in y_axis_value:
        if i < 100:
            noinf.append(i)

    x_threshold = 0.584963
    y_threshold = 1.5

    ymin = 0
    ymax = math.ceil(max(noinf))
    xmin = -4
    xmax = 4

    df['group'] = 'darkslategrey'
    df.loc[(x_axis_value > x_threshold) & (y_axis_value > y_threshold), 'group'] = 'maroon'
    df.loc[(x_axis_value < -x_threshold) & (y_axis_value > y_threshold), 'group'] = 'deepskyblue'
    df.loc[y_axis_value < y_threshold, 'group'] = 'lightgrey'
    df.columns = ['ipi', 'description'] + a + ['log2FC', 'p value', '-log10FDR', 'change']

    df.to_csv('{}_{}_combined_volcano.txt'.format(basename, query),
              index=True, header=True, sep='\t',
              )

    ax.vlines(-x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)
    ax.vlines(x_threshold, ymin, ymax, color='dimgrey', linestyle='dashed', linewidth=1)
    ax.hlines(y_threshold, xmin, xmax, color='dimgrey', linestyle='dashed', linewidth=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

    params = {'family':'sans-serif',
              'sans-serif': 'Helvetica',
              'weight': 'normal'}
    plt.rc('font', **params)

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel('${log_2(Fold\ change)}$', fontsize = 19)
    plt.ylabel('${-log_{10}(Q\ value)}$', fontsize = 19)
    plt.xticks(np.arange(xmin, xmax+1, step=1),fontsize=21)
    plt.yticks(np.arange(ymin, ymax+1, step=1),fontsize=21)
    plt.tick_params(axis='both', length = 10, width = 2)

    ax.scatter(x_axis_value, y_axis_value, s=23, c=df['change'], alpha=0.7)

    plt.savefig('{}_{}_volcano.pdf'.format(basename, query))
    #plt.show()

    endtime = datetime.datetime.now()
    print('time use is {} seconds'.format((endtime - starttime).seconds))

def combination_filter(basename, query):

    file_list = get_file(basename)
    query = query
    compiled_ms_info = {}
    common_keys = []
    i = 0
    while i < len(file_list):
        compiled_ms_info[i] = read_single_combined_dta(file_list[i])
        common_keys.append(compiled_ms_info[i].keys())
        i += 1

    # Finding common consensus sequences
    a = common_keys[0] & common_keys[1]
    i = 2
    while i < len(file_list):
        a = a & common_keys[i]
        i += 1

    i = 0
    combined_ms_info = defaultdict(list)
    # Adding UniprotID, gene_description
    for consensus_seq in a:
        combined_ms_info[consensus_seq] = [0,0]
        combined_ms_info[consensus_seq][0] = compiled_ms_info[0][consensus_seq][0]
        combined_ms_info[consensus_seq][1] = compiled_ms_info[0][consensus_seq][1]
    # UniprotID, gene_description, L/H #1, L/M #1, L/H #2, L/M #2, L/H #3, L/M #3,
    while i < len(file_list):
        for consensus_seq in a:
            combined_ms_info[consensus_seq] += compiled_ms_info[i][consensus_seq][2:4]
        i += 1

    # Calculate FDR
    #print(filtered_ms_info(combined_ms_info, query).items())
    sorted_fdr_dict = sorted(filtered_ms_info(combined_ms_info, query).items(), key=lambda x: x[-1][-1])
    #print(sorted_fdr_dict)
    len_dict = len(sorted_fdr_dict)
    fdred_dict = {}
    for i in sorted_fdr_dict:
        adjusted_q_value = -1 * np.log10 ((i[1][-1]*len_dict)/(sorted_fdr_dict.index(i) + 1))
        fdred_dict[i[0]] = i[1] + [adjusted_q_value]
    df = pd.DataFrame(pd.DataFrame(fdred_dict).T)
    return df

def filtered_ms_info(combined_ms_info, query):
    filtered_ms_infom = {}
    if query == 'heavy':
        for key in combined_ms_info.keys():
            n_list = combined_ms_info[key][0:2] + combined_ms_info[key][2::2]
            if retun_info_log2_p(n_list) != None:
                filtered_ms_infom[key] = retun_info_log2_p(n_list)
    else:
        for key in combined_ms_info.keys():
            n_list = combined_ms_info[key][0:2] + combined_ms_info[key][3::2]
            if retun_info_log2_p(n_list) != None:
                filtered_ms_infom[key] = retun_info_log2_p(n_list)
    #print(filtered_ms_infom)
    return filtered_ms_infom

def retun_info_log2_p(four_ele_list):
    # Remove UniprotID contains Reversed
    # Filter L/H in both reps which value is 0
    #print(four_ele_list)
    cal_value = four_ele_list[2::]
    # Exclude [15, 15, 15 ...] or [0.7, 0.7, 0.7 ...] list
    if np.prod(cal_value) != 0 and ('Rev' in four_ele_list[0]) == False \
            and (list(set(cal_value)) in [[15.0],[0.07]] ) == False:
#        print(four_ele_list)
        retun_info_log2_p_list = recognize_elements_legal(four_ele_list)
        return retun_info_log2_p_list

def recognize_elements_legal(four_ele_list):
    cal_value = four_ele_list[2::]
    if len(set(cal_value)) == 1:
        # this rel_cal_value is generated to avoid
        # [A, A], [1, 1] ==> [A-0.01, A], [1, 1]
        # test1: rel_cal_value = [cal_value[0] - 0.01] + cal_value[1::]
        rel_cal_value = cal_value + [cal_value[0] - 0.01]
        p_value, log2 = ret_p_value_log2(rel_cal_value)
        retun_info_log2_p_list = four_ele_list[0:2] + rel_cal_value[0:-1] + [log2] + [p_value]
    else:
        rel_cal_value = cal_value
        p_value, log2 = ret_p_value_log2(rel_cal_value)
        retun_info_log2_p_list = four_ele_list[0:2] + rel_cal_value + [log2] + [p_value]
    return retun_info_log2_p_list

def ret_p_value_log2(list_a):
    p_value = float(scipy.stats.ttest_ind(list_a, [1 for _ in range(len(list_a))])[1])  # to avoid [1,1],[1,1]
    log2 = float(np.log2(np.mean(list_a)))
    return p_value, log2

def read_single_combined_dta(file):
    file_input = open(file)
    ms_info = defaultdict(list)
    ipi_list = []
    descrip_list = []
    consensus_sequence = []
    light_heavy = []
    light_medium = []
    for line in file_input:
        line = line.strip().split('\t')
        if str.isdigit(line[0]) == True:
            ipi_list.append(0)
            descrip_list.append(0)
            consensus_sequence.append(line[4])
            light_heavy.append(float(line[6]))
            light_medium.append(float(line[7]))
        else:
            ipi_list.append(line[0])
            descrip_list.append(line[1])
    gen_list_ipi_list = gen_list(ipi_list)
    gen_list_descrip_list = gen_list(descrip_list)

    # Generation of a dictionary
    # Key is consensus sequence
    for i in consensus_sequence:
        index = consensus_sequence.index(i)
        ms_info[i] = [gen_list_ipi_list[index],
                      gen_list_descrip_list[index],
                      light_heavy[index],
                      light_medium[index]]
    return ms_info

def gen_list(ipi_list):
    tis = [ipi_list[i+1] for i, v in enumerate(ipi_list) if v == 0 ]
    return tis

########################################

def get_file(base_name):
    file_list = []
    current_work_dir = os.getcwd()
    all_file = os.listdir(current_work_dir)
    for each_file in all_file:
        if base_name in each_file and (('vol' in each_file) == False) and os.path.splitext(each_file)[-1] == '.txt':
            file_list.append(each_file)
    print(file_list)
    return file_list

plot_volcano(basename=args.basename,
             query=args.query
             )