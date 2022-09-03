from Bio import SeqIO
import re
import os
import RNA
import errors
import configparser
from subprocess import Popen
import pandas as pd
from difflib import SequenceMatcher
import time

# GAGGUAG - let7 seed
# CACCGGG - mir35 seed

# constant variables:
dirpath = os.getcwd()

# default parameters for the pipeline
hairpinFile = ''
matureFile = ''
# output_file_all = "output_all.csv"
output_file_name_filtered_by_seed = "output_filtered_by_seed.csv"
seed = 'CACCGGG'

# read settings file
settings = configparser.ConfigParser()
settings._interpolation = configparser.ExtendedInterpolation()
settings.read('settings.ini')

if not settings.has_section('mode_0'):
    print("mode 0 error: " + errors.errorMessage('section_not_exist'))
if settings.has_option('mode_0', 'hairpin_file_path'):
    hairpinFile = settings.get('mode_0', 'hairpin_file_path')
else:
    print("mode 0 error: " + errors.errorMessage('missing_input_file'))
if settings.has_option('mode_0', 'mature_file_path'):
    matureFile = settings.get('mode_0', 'mature_file_path')
else:
    print("mode 0 error: " + errors.errorMessage('missing_input_file'))

# if settings.has_option('mode0', 'output_file_name_all'):
#     output_file_all = settings.get('mode0', 'output_file_name_all')
if settings.has_option('mode_0', 'output_file_name_filtered_by_seed'):
    output_file_name_filtered_by_seed = settings.get('mode_0', 'output_file_name_filtered_by_seed')
if settings.has_option('mode_0', 'seed'):
    seed = settings.get('mode_0', 'seed')


# functions
def ct_file_parser_3p(ct_df, end_hairpin):
    index_i = end_hairpin
    repair_index = 0
    decreased = False
    while int(ct_df.loc[index_i][4]) == 0:
        index_i -= 1
        repair_index += 1
        decreased = True
    index_j = int(ct_df.loc[index_i][4])
    if decreased:
        return max(0, index_j - repair_index), end_hairpin
    return index_j, end_hairpin


def ct_file_parser_5p(ct_df, start_hairpin):
    index_i = start_hairpin
    repair_index = 0
    increased = False
    while int(ct_df.loc[index_i][4]) == 0:
        index_i += 1
        repair_index += 1
        increased = True
    index_j = int(ct_df.loc[index_i][4])
    if increased:
        return start_hairpin, min(len(ct_df), index_j + repair_index)
    return start_hairpin, index_j


def find_loop_size_3p(ct_df, start_mature):
    if int(ct_df.loc[start_mature - 1][4]) != 0:
        start_loop = int(ct_df.loc[start_mature - 1][4]) + 2  # Overhead
    else:
        index = start_mature + 1
        repair_index = 2
        while int(ct_df.loc[index][4]) == 0:
            index += 1
            repair_index += 1
        start_loop = int(ct_df.loc[index][4]) + repair_index

    end_loop = start_mature

    return start_loop, end_loop


def find_loop_size_5p(ct_df, end_mature):
    if int(ct_df.loc[end_mature + 1][4]) != 0:
        end_loop = int(ct_df.loc[end_mature + 1][4]) + 2  # Overhead
    else:
        index = end_mature - 1
        repair_index = 2
        while int(ct_df.loc[index][4]) == 0:
            index -= 1
            repair_index += 1
        end_loop = int(ct_df.loc[index][4]) - repair_index

    start_loop = end_mature

    return start_loop, end_loop


def create_csv_data():
    write_file_seed = open(output_file_name_filtered_by_seed, "w")
    try:
        i = 1
        j = 1
        res = {}
        # csv headlines:
        string_file_to_write_seed = "id,organism,hairpin_id,hairpin_seq,seq_cut_tails,cutted_hairpin_length,fold,energy,\
                                    mature_3p,mature_5p,numbers_of_connections,bp_ratio,loop_size,star,star_length" + '\n'

        records_dict_mature = list(SeqIO.parse(matureFile, "fasta"))
        records_dict_hairpin = list(SeqIO.parse(hairpinFile, "fasta"))

        # loop over all matures and find match/es in hairpin file
        for record_h in records_dict_hairpin:
            if j % 5000 == 0:
                print(j)
            j = j + 1

            # filter dict_hairpin by seed -> only hairpins with the given seed stay
            if seed not in str(record_h.seq):
                continue
            # if j > 100:
            #     break
            for record_m in records_dict_mature:
                # get only the right organism by 3 letters (for example: c.elegnas = cel)
                if record_m.name.split('-')[0] != record_h.name.split('-')[0]:
                    continue
                # find matched mature in hairpin
                m = re.search(str(record_m.seq), str(record_h.seq))
                if m == None:
                    continue
                # check if 3p or 5p
                if m.start() > (len(record_h.seq) - m.end()):
                    mature3p = str(record_m.seq)
                    mature5p = ""
                else:
                    mature5p = str(record_m.seq)
                    mature3p = ""

                # extract organism name from description
                organism = record_m.description.split(' ')[2] + " " + record_m.description.split(' ')[3]

                # get the fold from hairpin
                fold = RNA.fold(str(record_h.seq))
                # hairpin = str(record_h.seq)

                # add record to records dict
                if mature3p != '':
                    key = record_h.name + '-3p'
                    if key not in res:
                        res[key] = {}
                    else:
                        match_exist = SequenceMatcher(None, record_h.name, res[key]['mature_name'])
                        match_new = SequenceMatcher(None, record_h.name, record_m.name)
                        if match_exist.ratio() > match_new.ratio():
                            continue


                else:  # if mature5p != '':
                    key = record_h.name + '-5p'
                    if record_h.name not in res:
                        res[key] = {}
                    else:
                        match_exist = SequenceMatcher(None, record_h.name, res[key]['mature_name'])
                        match_new = SequenceMatcher(None, record_h.name, record_m.name)
                        if match_exist.ratio() > match_new.ratio():
                            continue

                res[key]['hairpin_name'] = record_h.name
                res[key]['mature_name'] = record_m.name
                res[key]['hairpin'] = str(record_h.seq)
                res[key]['fold'] = fold[0]
                res[key]['organism'] = organism
                res[key]['mature3p'] = mature3p
                res[key]['mature5p'] = mature5p


        write_file_seed.write(string_file_to_write_seed)
        for k, v in res.items():

            if seed in v['mature3p'] or seed in v['mature5p']:
                if seed in v['mature3p']:
                    if seed != v['mature3p'][1:8]:
                        continue
                else:
                    if seed != v['mature5p'][1:8]:
                        continue
                with open('ct.txt', 'w') as infile:
                    infile.write('>' + str(k) + '\n')
                    infile.write(v['hairpin'] + '\n')
                    infile.write(v['fold'])

                cmd = "RNAfold --noPS ct.txt | b2ct > ct_file.ct"
                p = Popen(cmd, shell=True)
                p.communicate()

                ct_df = pd.read_csv('ct_file.ct', delimiter='\s+', header=None, names=[0, 1, 2, 3, 4, 5])
                ct_df = ct_df.iloc[1:]
                ct_df.astype({4: 'int'}).dtypes

                # find indexes of the seed, mature, and hairpin
                start_seed, end_seed = find_seed(seed, v['hairpin'])

                # in 3p and in 5p is different?
                start_mature = start_seed - 1
                end_mature = min(end_seed + 14, len(ct_df))

                numbers_of_connections = mature_complimentarity(ct_df.loc[start_mature:end_mature])

                if seed in v['mature3p']:
                    bp_ratio = numbers_of_connections / float(len(v['mature3p']))

                    # find indexes of start and end of the hairpin from ct file
                    hairpin_boundries = ct_file_parser_3p(ct_df, end_mature)

                    # find the indexes of the loop
                    start_loop, end_loop = find_loop_size_3p(ct_df, start_mature)

                    start_star = hairpin_boundries[0] + 2
                    end_star = start_loop - 1
                    star_length = end_star - start_star

                    star = v['hairpin'][start_star - 1:end_star]


                else:
                    bp_ratio = numbers_of_connections / float(len(v['mature5p']))
                    # find indexes of start and end of the hairpin from ct file
                    hairpin_boundries = ct_file_parser_5p(ct_df, start_mature)
                    # find the indexes of the loop
                    start_loop, end_loop = find_loop_size_5p(ct_df, end_mature)

                    start_star = hairpin_boundries[1] - 2
                    end_star = end_loop + 1
                    star_length = start_star - end_star

                    star = v['hairpin'][end_star + 1:start_star]

                if hairpin_boundries[0] > hairpin_boundries[1]:
                    continue

                # cut the hairpin with the new indexes
                cutted_hairpin = v['hairpin'][hairpin_boundries[0]:hairpin_boundries[1]]

                loop_size = end_loop - start_loop
                fold_for_energy = RNA.fold(cutted_hairpin)

                write_file_seed.write(
                    str(i) + "," + v['organism'] + "," + k + "," + v['hairpin'] + "," + cutted_hairpin + "," +
                    str(len(cutted_hairpin)) + "," + v['fold'] + "," + str(fold_for_energy[1]) + "," + v['mature3p'] +
                    "," + v['mature5p'] + "," + str(numbers_of_connections) + "," + str(bp_ratio) + "," +
                    str(loop_size) + "," + star + "," + str(star_length) + '\n')

            i = i + 1
        write_file_seed.close()

    finally:
        if not write_file_seed.closed:
            write_file_seed.close()


def find_seed(seed, seq):
    for m in re.finditer(seed, seq):
        start = m.start() + 1
        end = m.end()
        return start, end


def mature_complimentarity(mature_df):
    mature_connections = 0
    for index, row in mature_df.iterrows():
        if int(row[4]) != 0:
            # check only connection outside mature
            if mature_df[0].iloc[0] <= int(row[4]) <= mature_df[0].iloc[len(mature_df)-1]:
                continue
            mature_connections += 1
    return mature_connections


# add energy to candidates by final fold:
def get_fold_energy_for_candidate(fold_seq):
    if len(fold_seq) > 0:
        fold = RNA.fold(fold_seq)
        return str(fold[1])
    return ''


def reverse_hairpin(hairpin_fold):
    hairpin_fold_r = hairpin_fold[::-1]
    start_replace_from = hairpin_fold_r.count(')')
    hairpin_fold_r = hairpin_fold_r.replace("(", ")")
    hairpin_fold_r = hairpin_fold_r.replace(")", "(", start_replace_from)
    return hairpin_fold_r





print("start running mode 0..")
start = time.clock()
create_csv_data()
elapsed = (time.clock() - start)
print("Program executed in " + str(elapsed))
print("finish running mode 0")





#
# records_dict_mature = list(SeqIO.parse(matureFile, "fasta"))
# records_dict_hairpin = list(SeqIO.parse(hairpinFile, "fasta"))
# # loop over all matures and find match/es in hairpin file
# for record_m in records_dict_mature:
#     for record_h in records_dict_hairpin:
#         m = re.search(str(record_m.seq), str(record_h.seq))
#         if m == None:
#             continue
#         if record_h.name == 'cel-mir-38':
#             print("yey")