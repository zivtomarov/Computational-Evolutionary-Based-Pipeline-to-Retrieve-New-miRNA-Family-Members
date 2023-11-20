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
import subprocess
from io import StringIO


# GAGGUAG - let7 seed
# CACCGGG - mir35 seed

# constant variables:
dirpath = os.getcwd()

# default parameters for the pipeline
hairpinFile = ''
matureFile = ''
# output_file_all = "output_all.csv"
output_file_name_filtered_by_seed = "output_filtered_by_seed.csv"
# seed = 'CACCGGG'

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

# Nematodes from mirBase - the parameters from positive DB regarding these families only
if settings.has_option('mode_0', 'organism_list'):
    organism_list = settings['mode_0']['organism_list'].split(',')


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

def ct_file_parser_3p_check(ct_df, start_mature, end_mature, param=30):
    index_i = end_mature
    repair_index_end_mature = 0
    decreased_end_mature = False
    valid = True
    while int(ct_df.loc[index_i][4]) == 0:
        index_i -= 1
        repair_index_end_mature += 1
        decreased_end_mature = True
    start_hairpin = int(ct_df.loc[index_i][4])-1  # file starts from index 1, so  decrease 1 for index to align with index 0
    repair_index_start_star = repair_index_end_mature


    index_i = start_mature
    repair_index_start_mature = 0
    decreased_start_mature = False
    while int(ct_df.loc[index_i][4]) == 0:
        index_i += 1
        repair_index_start_mature += 1
        decreased_start_mature = True
    # end_star = int(ct_df.loc[index_i][4])-1  # file starts from index 1, so  decrease 1 for index to align with index 0

    direct = False
    if int(ct_df.loc[start_mature - 2][4]) != 0:
        end_star_direct = int(ct_df.loc[start_mature - 2][4]) - 1
        direct = True
    else:
        end_star_undirect = int(
            ct_df.loc[index_i][4]) - 1  # file starts from index 1, so  decrease 1 for index to align with index 0

    if start_hairpin > end_mature or start_hairpin > param:
        valid = False

        if decreased_start_mature:
            if direct:
                return {'start_hairpin': max(0, start_hairpin - repair_index_end_mature) + 2,
                        'end_hairpin': end_mature - 1,
                        'start_star': max(0, start_hairpin - repair_index_start_star) + 2,
                        'end_star': end_star_direct, 'valid': valid}
            return {'start_hairpin': max(0, start_hairpin - repair_index_end_mature) + 2, 'end_hairpin': end_mature - 1,
                    'start_star': max(0, start_hairpin - repair_index_start_star) + 2,
                    'end_star': end_star_undirect + repair_index_start_mature + 2, 'valid': valid}
        if direct:
            return {'start_hairpin': max(0, start_hairpin - repair_index_end_mature) + 2, 'end_hairpin': end_mature - 1,
                    'start_star': max(0, start_hairpin - repair_index_start_star + 2),
                    'end_star': end_star_direct, 'valid': valid}
        return {'start_hairpin': max(0, start_hairpin - repair_index_end_mature) + 2, 'end_hairpin': end_mature - 1,
                'start_star': max(0, start_hairpin - repair_index_start_star + 2),
                'end_star': end_star_undirect + 2, 'valid': valid}
    if decreased_start_mature:
        if direct:
            return {'start_hairpin': start_hairpin + 2, 'end_hairpin': end_mature - 1,
                    'start_star': max(0, start_hairpin - repair_index_start_star) + 2,
                    'end_star': end_star_direct, 'valid': valid}
        return {'start_hairpin': start_hairpin + 2, 'end_hairpin': end_mature - 1,
                'start_star': max(0, start_hairpin - repair_index_start_star) + 2,
                'end_star': end_star_undirect + 2, 'valid': valid}
    if direct:
        return {'start_hairpin': start_hairpin + 2, 'end_hairpin': end_mature - 1, 'start_star': start_hairpin + 2,
                'end_star': end_star_direct, 'valid': valid}
    return {'start_hairpin': start_hairpin, 'end_hairpin': end_mature - 1, 'start_star': start_hairpin + 2,
            'end_star': end_star_undirect + 2, 'valid': valid}


def ct_file_parser_5p_check(ct_df, start_mature, end_mature, param):
    index_i = start_mature
    repair_index_start_mature = 0
    increased_start_mature = False
    valid = True
    while int(ct_df.loc[index_i][4]) == 0:
        index_i += 1
        repair_index_start_mature += 1
        increased_start_mature = True
    end_hairpin = int(ct_df.loc[index_i][4])-1
    repair_index_end_star = repair_index_start_mature

    index_i = end_mature
    repair_index_end_mature = 0
    increased_end_mature = False
    while int(ct_df.loc[index_i][4]) == 0:
        index_i -= 1
        repair_index_end_mature += 1
        increased_end_mature = True
    # start_star = int(ct_df.loc[index_i][4])-1

    direct = False
    if int(ct_df.loc[end_mature - 2][4]) != 0:
        start_star_direct = int(ct_df.loc[end_mature - 2][4]) - 1
        direct = True
    else:
        start_star_undirect = int(
            ct_df.loc[index_i][4]) - 1  # file starts from index 1, so  decrease 1 for index to align with index 0

    # fix problem that end of star at 5p is after the calculated end_hairpin
    if min(len(ct_df), end_hairpin + repair_index_start_mature) < min(len(ct_df), end_hairpin + repair_index_end_star+2):
        repair_index_start_mature = repair_index_end_star+2

    if end_hairpin < start_mature or end_hairpin < param:
        valid = False

    if increased_start_mature:
        if increased_end_mature:
            if direct:
                return {'start_hairpin': start_mature - 1,
                        'end_hairpin': min(len(ct_df), end_hairpin + repair_index_start_mature) + 2,
                        'start_star': start_star_direct,
                        'end_star': min(len(ct_df), end_hairpin + repair_index_end_star) + 2, 'valid': valid}
            return {'start_hairpin': start_mature - 1,
                    'end_hairpin': min(len(ct_df), end_hairpin + repair_index_start_mature),
                    'start_star': start_star_undirect - repair_index_end_mature + 2,
                    'end_star': min(len(ct_df), end_hairpin + repair_index_end_star) + 2, 'valid': valid}
        if direct:
            return {'start_hairpin': start_mature - 1,
                    'end_hairpin': min(len(ct_df), end_hairpin + repair_index_start_mature) + 2,
                    'start_star': start_star_direct,
                    'end_star': min(len(ct_df), end_hairpin + repair_index_end_star) + 2, 'valid': valid}
        return {'start_hairpin': start_mature - 1,
                'end_hairpin': min(len(ct_df), end_hairpin + repair_index_start_mature) + 2,
                'start_star': start_star_undirect + 2,
                'end_star': min(len(ct_df), end_hairpin + repair_index_end_star) + 2, 'valid': valid}
    if increased_end_mature:
        if direct:
            return {'start_hairpin': start_mature - 1,
                    'end_hairpin': min(len(ct_df), end_hairpin + repair_index_start_mature) + 2,
                    'start_star': start_star_direct, 'end_star': end_hairpin + 2, 'valid': valid}
        return {'start_hairpin': start_mature - 1,
                'end_hairpin': min(len(ct_df), end_hairpin + repair_index_start_mature) + 2,
                'start_star': start_star_undirect - repair_index_end_mature + 2, 'end_star': end_hairpin + 2,
                'valid': valid}
    if direct:
        return {'start_hairpin': start_mature - 1, 'end_hairpin': end_hairpin, 'start_star': start_star_direct,
                'end_star': end_hairpin + 2, 'valid': valid}

    return {'start_hairpin': start_mature - 1, 'end_hairpin': end_hairpin, 'start_star': start_star_undirect + 2,
            'end_star': end_hairpin + 2, 'valid': valid}

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
        # string_file_to_write_seed = "ID,Organism,hairpin_name,Hairpin_seq,Hairpin_seq_trimmed,Hairpin_seq_trimmed_length,Fold," \
        #                             "Energy,mature_3p,mature_5p,numbers_of_connections,bp_ratio,loop_size,star," \
        #                             "star_length" + '\n'
        string_file_to_write_seed = "ID,Organism,Hairpin_name,Hairpin_seq,Hairpin_seq_trimmed," \
                                    "Hairpin_seq_trimmed_length,Fold," \
                                    "Energy,Mature,Mature_connections,Mature_BP_ratio,Mature_max_bulge,Loop_length,Star," \
                                    "Star_length,Star_connections,Star_BP_ratio,Star_max_bulge,Max_bulge_symmetry,Max_one_mer_mature,"\
                                    "Max_two_mer_mature,Max_one_mer_hairpin,Max_two_mer_hairpin,5p/3p" + '\n'

        records_dict_mature = list(SeqIO.parse(matureFile, "fasta"))
        records_dict_hairpin = list(SeqIO.parse(hairpinFile, "fasta"))

        # loop over all matures and find match/es in hairpin file
        for record_h in records_dict_hairpin:
            # if j % 100 == 0:
            #     break
            if j % 5000 == 0:
                print(j)
            j = j + 1

            # filter dict_hairpin by seed -> only hairpins with the given seed stay
            if seed not in str(record_h.seq):
                continue
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
            mature_3p = False
            mature_5p = False
            if seed in v['mature3p'] or seed in v['mature5p']:
                if seed in v['mature3p']:
                    if seed != v['mature3p'][1:8]:
                        continue
                else:
                    if seed != v['mature5p'][1:8]:
                        continue
                # with open('ct.txt', 'w') as infile:
                #     infile.write('>' + str(k) + '\n')
                #     infile.write(v['hairpin'] + '\n')
                #     infile.write(v['fold'])
                #
                # cmd = "RNAfold --noPS ct.txt | b2ct > ct_file.ct"
                # p = Popen(cmd, shell=True)
                # p.communicate()

                ct_fold = v['fold']
                ct_hairpin = v['hairpin']
                ct_data = f'>{key}\n{ct_hairpin}\n{ct_fold}'
                cmd = "RNAfold --noPS | b2ct"
                p = subprocess.run(cmd, shell=True, input=ct_data, universal_newlines=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
                out, err = p.stdout, p.stderr


                # ct_df = pd.read_csv('ct_file.ct', delimiter='\s+', header=None, names=[0, 1, 2, 3, 4, 5])
                ct_df = pd.read_csv(StringIO(out), delimiter='\s+', header=None, names=[0, 1, 2, 3, 4, 5])[1:]
                ct_df = ct_df.astype({4: 'int'})

                # ct_df = ct_df.iloc[1:]
                # ct_df.astype({4: 'int'}).dtypes

                # find indexes of the seed, mature, and hairpin
                start_seed, end_seed = find_seed(seed, v['hairpin'])


                start_mature = start_seed - 1
                end_mature = min(end_seed + 14, len(ct_df))

                mature_df = ct_df.loc[start_mature:end_mature]
                mature_numbers_of_connections, mature_max_bulge = mature_complimentarity(mature_df)

                mature_bp_ratio = mature_numbers_of_connections / float(len(mature_df))

                ## deleted to check the parsers

                # if seed in v['mature3p']:
                #     # bp_ratio = numbers_of_connections / float(len(v['mature3p']))
                #
                #     # find indexes of start and end of the hairpin from ct file
                #     hairpin_boundries = ct_file_parser_3p(ct_df, end_mature)
                #
                #     # find the indexes of the loop
                #     start_loop, end_loop = find_loop_size_3p(ct_df, start_mature)
                #
                #     start_star = hairpin_boundries[0] + 2
                #     end_star = start_loop - 1
                #     star_length = end_star - start_star
                #
                #     star = v['hairpin'][start_star - 1:end_star]
                #
                #     mature_3p = True
                #
                #
                # else:
                #     # bp_ratio = numbers_of_connections / float(len(v['mature5p']))
                #
                #     # find indexes of start and end of the hairpin from ct file
                #     hairpin_boundries = ct_file_parser_5p(ct_df, start_mature)
                #
                #     # find the indexes of the loop
                #     start_loop, end_loop = find_loop_size_5p(ct_df, end_mature)
                #
                #     start_star = hairpin_boundries[1] - 2
                #     end_star = end_loop + 1
                #     star_length = start_star - end_star
                #
                #     star = v['hairpin'][end_star + 1:start_star]
                #
                #     mature_5p = True

                ###############################################################
                if seed in v['mature3p']:
                    hairpin_boundries = ct_file_parser_3p_check(ct_df, start_mature, end_mature, 30)
                    # start_star = hairpin_boundries['start_star']
                    # end_star = hairpin_boundries['end_star']
                    # star_length = end_star - start_star + 1
                    # star = v['hairpin'][start_star:end_star + 1]

                    start_loop, end_loop = find_loop_size_3p(ct_df,
                                                             start_mature)

                    mature_3p = True

                else:
                    hairpin_boundries = ct_file_parser_5p_check(ct_df, start_mature, end_mature, 30)
                    # start_star = hairpin_boundries['start_star']
                    # end_star = hairpin_boundries['end_star']
                    # star_length = end_star - start_star + 1
                    # star = v['hairpin'][start_star:end_star + 1]

                    start_loop, end_loop = find_loop_size_5p(ct_df,
                                                             end_mature)

                    mature_5p = True

                if hairpin_boundries['valid'] is False:
                    # del res[key]
                    # print("not valid")
                    continue

                start_star = hairpin_boundries['start_star']
                end_star = hairpin_boundries['end_star']
                star_length = end_star - start_star + 1
                star = v['hairpin'][start_star:end_star + 1]

                star_df = ct_df.loc[start_star+1:end_star+1]
                star_numbers_of_connections, star_max_bulge = star_complimentarity(star_df)

                star_bp_ratio = star_numbers_of_connections / float(len(star_df))

                max_bulge_symmetry = find_max_bulge_symmetry(mature_df)


                # cut the hairpin with the new indexes
                cutted_hairpin = v['hairpin'][hairpin_boundries['start_hairpin']:hairpin_boundries['end_hairpin'] + 1]

                ###############################################################

                ## deleted to check the parsers
                # if hairpin_boundries[0] > hairpin_boundries[1]:
                #     continue

                ## deleted to check the parsers
                # cut the hairpin with the new indexes
                # cutted_hairpin = v['hairpin'][hairpin_boundries[0]:hairpin_boundries[1]]

                mature = mature_df[1].str.cat()

                loop_size = end_loop - start_loop
                fold_and_energy = RNA.fold(cutted_hairpin)

                max_one_mer_mature = count_kmers(mature, 1) / len(mature)
                max_two_mer_mature = count_kmers(mature, 2) / len(mature)
                max_one_mer_hairpin = count_kmers(cutted_hairpin, 1) / len(cutted_hairpin)
                max_two_mer_hairpin = count_kmers(cutted_hairpin, 2) / len(cutted_hairpin)



                if mature_3p:
                    write_file_seed.write(
                        str(i) + "," + v['organism'] + "," + k + "," + v['hairpin'] + "," + cutted_hairpin + "," +
                        str(len(cutted_hairpin)) + "," + str(fold_and_energy[0]) + "," + str(fold_and_energy[1]) + "," +
                        mature + "," + str(mature_numbers_of_connections) + "," + str(mature_bp_ratio) + "," +
                        str(mature_max_bulge) + "," + str(loop_size) + "," + star + "," + str(star_length) + "," +
                        str(star_numbers_of_connections) + "," + str(star_bp_ratio) + "," + str(star_max_bulge) + "," +
                        str(max_bulge_symmetry) + "," + str(max_one_mer_mature) + "," + str(max_two_mer_mature) + "," +
                        str(max_one_mer_hairpin) + "," + str(max_two_mer_hairpin) + "," + "3p"'\n')
                elif mature_5p:
                    write_file_seed.write(
                        str(i) + "," + v['organism'] + "," + k + "," + v['hairpin'] + "," + cutted_hairpin + "," +
                        str(len(cutted_hairpin)) + "," + str(fold_and_energy[0]) + "," + str(fold_and_energy[1]) + "," +
                        mature + "," + str(mature_numbers_of_connections) + "," + str(mature_bp_ratio) + "," +
                        str(mature_max_bulge) + "," + str(loop_size) + "," + star + "," + str(star_length) + "," +
                        str(star_numbers_of_connections) + "," + str(star_bp_ratio) + "," + str(star_max_bulge) + "," +
                        str(max_bulge_symmetry) + "," + str(max_one_mer_mature) + "," + str(max_two_mer_mature) + "," +
                        str(max_one_mer_hairpin) + "," + str(max_two_mer_hairpin) + "," +"5p"'\n')

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


def find_max_bulge_symmetry(mature_df):

    mature_i = 0
    star_i = 0
    max_bulge_symmetry = 0
    seen_bulge = False

    # skip all first rows with zero connections - not a bulge
    counter = 0
    for row in mature_df.values:
        if int(row[4]) == 0:
            counter += 1
            continue
        else:
            break

    for index, row in enumerate(mature_df.values[counter:]):
        # if index < counter:
        #     continue
        if int(row[4]) == 0 and not seen_bulge:
            seen_bulge = True
            mature_i = row[0] - 1  # start of bulge
        else:
            if int(row[4]) == 0:
                continue

            elif not seen_bulge:
                star_i = int(row[4])
                continue

            else:
                seen_bulge = False

                mature_j = row[0]
                star_j = int(row[4])

                mature_bulge = abs(mature_j - mature_i) - 1
                star_bulge = abs(star_j - star_i) - 1
                diff_bulge = abs(mature_bulge - star_bulge)
                if diff_bulge > max_bulge_symmetry:
                    max_bulge_symmetry = diff_bulge

                star_i = star_j

    return max_bulge_symmetry


def mature_complimentarity(mature_df):
    bulge_flag = False
    max_bulge = 0
    count_bulge = 0
    mature_connections = 0
    for row in mature_df.values:
        if int(row[4]) != 0:
            if mature_df[0].iloc[0] <= int(row[4]) <= mature_df[0].iloc[len(mature_df) - 1]:  # check only connection outside mature
                continue
            mature_connections += 1

            if bulge_flag:
                max_bulge = max(max_bulge, count_bulge)
                # if max_bulge < count_bulge:
                #     max_bulge = count_bulge
                count_bulge = 0
                bulge_flag = False
        else:
            count_bulge += 1
            bulge_flag = True

    return mature_connections, max_bulge




def star_complimentarity(star_df):
    bulge_flag = False
    max_bulge = 0
    count_bulge = 0
    star_connections = 0
    for row in star_df.values:
        if row[4] != 0:
            # check only connection outside star
            if star_df[0].iloc[0] <= int(row[4]) <= star_df[0].iloc[len(star_df) - 1]:
                continue
            star_connections += 1

            if bulge_flag:
                max_bulge = max(max_bulge, count_bulge)
                # if max_bulge < count_bulge:
                #     max_bulge = count_bulge
                count_bulge = 0
                bulge_flag = False
        else:
            count_bulge += 1
            bulge_flag = True

    return star_connections, max_bulge


def count_kmers(sequence, k):
    """Count kmer occurrences in a given sequence.

    Parameters
    ----------
    read : string
        A single DNA sequence.
    k : int
        The value of k for which to count kmers.

    Returns
    -------
    counts : dictionary, {'string': int}
        A dictionary of counts keyed by their individual kmers (strings
        of length k).

    """
    # Start with an empty dictionary
    counts = {}
    # Calculate how many kmers of length k there are
    num_kmers = len(sequence) - k + 1
    # Loop over the kmer start positions
    for i in range(num_kmers):
        # Slice the string to get the kmer
        kmer = sequence[i:i+k]
        # Add the kmer to the dictionary if it's not there
        if kmer not in counts:
            counts[kmer] = 0
        # Increment the count for this kmer
        counts[kmer] += 1
    # Return the final counts
    return max(counts.values())



def get_optimal_parameters(organism_list=None, positive_db_csv=None):
    if settings.has_option('mode_0', 'output_filter_parameters'):
        output_filter_parameters = settings['mode_0']['output_filter_parameters']

    if positive_db_csv is None:
        if settings.has_option('mode_0', 'output_file_name_filtered_by_seed'):
            positive_db_csv = settings.get('mode_0', 'output_file_name_filtered_by_seed')

    df = pd.read_csv(positive_db_csv)

    if organism_list is not None:
        df = df.loc[df['Organism'].isin(organism_list)]
    print(len(df))
    max_energy = float("{:.2f}".format(df["Energy"].max()))
    min_energy = float("{:.2f}".format(df["Energy"].min()))
    min_star_length = float("{:.2f}".format(df["Star_length"].min()))
    max_star_length = float("{:.2f}".format(df["Star_length"].max()))
    min_loop_length = float("{:.2f}".format(df["Loop_length"].min()))
    max_loop_length = float("{:.2f}".format(df["Loop_length"].max()))
    min_mature_bp_ratio = float("{:.2f}".format(df["Mature_BP_ratio"].min()))
    max_mature_bp_ratio = float("{:.2f}".format(df["Mature_BP_ratio"].max()))
    min_star_bp_ratio = float("{:.2f}".format(df["Star_BP_ratio"].min()))
    max_star_bp_ratio = float("{:.2f}".format(df["Star_BP_ratio"].max()))
    min_trimmed_hairpin_length = float("{:.2f}".format(df["Hairpin_seq_trimmed_length"].min()))

    max_trimmed_hairpin_length = float("{:.2f}".format(df["Hairpin_seq_trimmed_length"].max()))
    max_mature_bulge = float("{:.2f}".format(df["Mature_max_bulge"].max()))
    max_star_bulge = float("{:.2f}".format(df["Star_max_bulge"].max()))
    max_bulge_symmetry = float("{:.2f}".format(df["Max_bulge_symmetry"].max()))

    max_one_mer_mature = float("{:.2f}".format(df["Max_one_mer_mature"].max()))
    max_two_mer_mature = float("{:.2f}".format(df["Max_two_mer_mature"].max()))
    max_one_mer_hairpin = float("{:.2f}".format(df["Max_one_mer_hairpin"].max()))
    max_two_mer_hairpin = float("{:.2f}".format(df["Max_two_mer_hairpin"].max()))

    dict_params = {'min_energy': min_energy, 'max_energy': max_energy,  'min_star_length': min_star_length,
                   'max_star_length': max_star_length, 'min_loop_length': min_loop_length,
                   'max_loop_length': max_loop_length, 'min_mature_bp_ratio': min_mature_bp_ratio,
                   'max_mature_bp_ratio': max_mature_bp_ratio, 'min_star_bp_ratio': min_star_bp_ratio,
                   'max_star_bp_ratio': max_star_bp_ratio, 'min_trimmed_hairpin_length': min_trimmed_hairpin_length,
                   'max_trimmed_hairpin_length': max_trimmed_hairpin_length, 'max_mature_bulge': max_mature_bulge,
                   'max_star_bulge': max_star_bulge, "max_bulge_symmetry": max_bulge_symmetry,
                   'max_one_mer_mature': max_one_mer_mature, 'max_two_mer_mature': max_two_mer_mature,
                   'max_one_mer_hairpin': max_one_mer_hairpin, 'max_two_mer_hairpin': max_two_mer_hairpin}

    f = open(output_filter_parameters, 'w')
    f.write(str(dict_params))
    f.close()

    # check and write to csv the ranges for the filter features for each organism from the nematodes list
    filter_features_ranges = df.groupby('Organism').agg({'Energy': ['mean', 'min', 'max'],
                                                         'Star_length': ['mean', 'min', 'max'],
                                                         'Loop_length': ['mean', 'min', 'max'],
                                                         'Mature_BP_ratio': ['mean', 'min', 'max'],
                                                         'Star_BP_ratio': ['mean', 'min', 'max'],
                                                         'Hairpin_seq_trimmed_length': ['mean', 'min', 'max'],
                                                         'Mature_max_bulge': ['mean', 'min', 'max'],
                                                         'Star_max_bulge': ['mean', 'min', 'max'],
                                                         'Max_bulge_symmetry': ['mean', 'min', 'max'],
                                                         'Max_one_mer_mature': ['mean', 'min', 'max'],
                                                         'Max_two_mer_mature': ['mean', 'min', 'max'],
                                                         'Max_one_mer_hairpin': ['mean', 'min', 'max'],
                                                         'Max_two_mer_hairpin': ['mean', 'min', 'max']})

    filter_features_ranges.to_csv("filter_features_ranges.csv")

    return dict_params


# Nematodes from mirBase - the parameters from positive DB regarding these families only
# organism_list = ['Ascaris suum', 'Brugia malayi', 'Caenorhabditis brenneri', 'Caenorhabditis briggsae',
#                  'Caenorhabditis elegans', 'Caenorhabditis remanei', 'Haemonchus contortus',
#                  'Heligmosomoides polygyrus', 'Pristionchus pacificus', 'Panagrellus redivivus', 'Strongyloides ratti']
# v = get_optimal_parameters(organism_list, 'output_let7_mode0_11_03_2022.csv')
# print("nematodes only: ", v)
#
# v = get_optimal_parameters(None, 'output_mir35_mode0_11_03_2022.csv')
# print("all positive DB: ", v)


# # add energy to candidates by final fold:
# def get_fold_energy_for_candidate(fold_seq):
#     if len(fold_seq) > 0:
#         fold = RNA.fold(fold_seq)
#         return str(fold[1])
#     return ''
#
#
# def reverse_hairpin(hairpin_fold):
#     hairpin_fold_r = hairpin_fold[::-1]
#     start_replace_from = hairpin_fold_r.count(')')
#     hairpin_fold_r = hairpin_fold_r.replace("(", ")")
#     hairpin_fold_r = hairpin_fold_r.replace(")", "(", start_replace_from)
#     return hairpin_fold_r


## TESTS
print("start running mode 0..")
start = time.clock()
# create_csv_data()
# v = get_optimal_parameters(organism_list, output_file_name_filtered_by_seed)
v = get_optimal_parameters(None, output_file_name_filtered_by_seed)

print("nematodes only parameters from positive DB: ", v)
elapsed = (time.clock() - start)
print("Program executed in " + str(elapsed))
print("finish running mode 0")
