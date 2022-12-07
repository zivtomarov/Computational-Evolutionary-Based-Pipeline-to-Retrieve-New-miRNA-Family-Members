import sys
import os
import time

import mode_1
import mode_2
import configparser
import errors
import pandas as pd


organism = []

# read settings file
settings = configparser.ConfigParser()
settings._interpolation = configparser.ExtendedInterpolation()
settings.read('settings.ini')


def checkIfExist(mode, organism1, organism2):

    if mode == 1:
        if settings.has_option('mode_1', 'output_path'):
            output_path = settings.get('mode_1', 'output_path')

            for filename in os.listdir(output_path):
                name = filename.split('.')[0]
                DBname = name.replace('results_', '')
                if organism1 == DBname:
                    print(f'{name} exists in folder')
                    return True
            print(f'{organism1} doesn''t exists in folder')
            return False

    elif mode == 2:
        if settings.has_option('mode_2', 'output_path'):
            output_path = settings.get('mode_2', 'output_path')
            for filename in os.listdir(output_path):
                name1 = filename.split(',')[0]
                name2 = filename.split(',')[1].split('.')[0]

                if (organism1 == name1 and organism2 == name2) or (organism1 == name2 and organism2 == name1):
                    return True
            return False

	# elif mode == 4:
	# 	if settings.has_option('mode4', 'output_path'):
	# 		output_path = settings.get('mode4', 'output_path')
    #
	# 		if os.path.exists(output_path + "/scores_mature_" + organism1):
	# 			for filename in os.listdir(output_path + "/scores_mature_" + organism1):
	# 				name1 = filename.split(',')[0].replace('scores_mature_', '')
	# 				name2 = filename.split(',')[1]
	# 				if (organism1 == name1 and organism2 == name2) or (organism1 == name2 and organism2 == name1):
	# 					return True
	# 				return False
    #
	# 		if os.path.exists(output_path + "/scores_mature_" + organism2):
	# 			for filename in os.listdir(output_path + "/scores_mature_" + organism2):
	# 				name1 = filename.split(',')[0].replace('scores_mature_', '')
	# 				name2 = filename.split(',')[1]
	# 				if (organism1 == name1 and organism2 == name2) or (organism1 == name2 and organism2 == name1):
	# 					return True
	# 				return False


def generate_mode1(organism_name_in_db):

    print(organism_name_in_db)
    start = time.clock()
    mode_1.find_candidates_by_seed()
    print("seed : ", len(mode_1.res))

    mode_1.filter_candidates()
    print("filter : ", len(mode_1.res))

    mode_1.add_fold_energy_to_candidates()
    mode_1.filterResultsByMaxEnergy()
    print("energy :", len(mode_1.res))

    mode_1.add_found_in_DB_to_res()
    print("TP :", len(mode_1.res))

    print("writing final results files")
    mode_1.write_final_results_file()
    mode_1.create_html_file()

    elapsed = (time.clock() - start)

    # clean the hashmap from mode 1 for the next candidate if exists
    mode_1.res = {}
    print("Program executed in " + str(elapsed))
    print("finish running mode 1")
    print("-----------------------------------------------------")

def generate_mode2():
    print("start running mode 2..")
    start = time.clock()
    try:
        if settings.has_option('mode_2', 'file_1_path'):
            file1_path = settings.get('mode_2', 'file_1_path')
        if settings.has_option('mode_2', 'file_2_path'):
            file2_path = settings.get('mode_2', 'file_2_path')
        if settings.has_option('mode_2', 'organism_file_1'):
            organism_file1 = settings.get('mode_2', 'organism_file_1')
        if settings.has_option('mode_2', 'organism_file_2'):
            organism_file2 = settings.get('mode_2', 'organism_file_2')

        if settings.has_option('mode_2', 'percentile_filter_percentage'):
            percentile_filter_percentage = float(settings.get('mode_2', 'percentile_filter_percentage'))
        if settings.has_option('mode_2', 'value_filter'):
            value_filter = int(settings.get('mode_2', 'value_filter'))
        if settings.has_option('mode_2', 'elbow_point'):
            elbow_point = settings.get('mode_2', 'elbow_point')

        # mode_2.align_files(file1_path, file2_path, organism_file1, organism_file2)
        mode_2.align_files(file1_path, file2_path, organism_file1, organism_file2, percentile_filter_percentage, value_filter, elbow_point)

        # modify indexes columns in df
		# arrange indexes correctly
        # saveIndexes(organism_file1, organism_file2)

        mode_2.find_number_of_orgs(file1_path, file2_path, organism_file1, organism_file2)

        elapsed = (time.clock() - start)
        print("Program executed in " + str(elapsed))
        print("finish running mode 2")
    except:
        print(errors.errorMessage('general_error'))
        print(sys.exc_info()[1])

def pipe():
    # # read settings file
    # settings = configparser.ConfigParser()
    # settings._interpolation = configparser.ExtendedInterpolation()
    # settings.read('settings.ini')

    if settings.has_option('mode_4', 'dir_fasta_path'):
        dir_fasta_path = settings.get('mode_4', 'dir_fasta_path')

        # go over all organisms fasta files to generate mode_1
        for filename in os.listdir(dir_fasta_path):
            name = filename.split('.')[0]
            DBname = name.replace('_', ' ')

            print("dir_fasta_path : " + dir_fasta_path)
            print("name : " + name)
            #
            # if name != 'caenorhabditis_remanei' and name != 'caenorhabditis_brenneri' and name != 'caenorhabditis_nigoni':
            #     continue

            if checkIfExist(1, name, 0):
                organism.append(name)
                continue

            # set the fasta file path to the next fasta file
            if settings.has_option('mode_1', 'fasta_file_path'):
                settings.set('mode_1', 'fasta_file_path', dir_fasta_path + '/' + filename)

            # set the output file name according to organism
            if settings.has_option('mode_1', 'output_file_name'):
                settings.set('mode_1', 'output_file_name', 'results_' + name)
                organism.append(name)

            # set organism name according to organism
            if settings.has_option('mode_1', 'organism_name_in_db'):
                settings.set('mode_1', 'organism_name_in_db', DBname[0].upper() + DBname[1:])

            # save changes to settings ini file
            with open('settings.ini', 'w') as settingFile:
                settings.write(settingFile)

            # Check that settings has been changed
            fasta_file_path = settings.get('mode_1', 'fasta_file_path')
            print("fasta_file_path : " + fasta_file_path)
            output_file_name = settings.get('mode_1', 'output_file_name')
            print("output_file_name : " + output_file_name)
            organism_name_in_db = settings.get('mode_1', 'organism_name_in_db')
            print("organism_name_in_db : " + organism_name_in_db)

            # Call mode 1 with new info
            # generate_mode1(organism_name_in_db)

        if settings.has_option('mode_1', 'output_path'):
            output_mode1_path = settings.get('mode_1', 'output_path')
        if settings.has_option('mode_4', 'dir_result_mode1_path'):
            settings.set('mode_4', 'dir_result_mode1_path', output_mode1_path)

        # Go over all possible couples of organisms from result files to generate MODE_2
        i = 0
        while i < len(organism):
            j = i
            while j < len(organism):
                # if (organism[i] == 'heligmosomoides_polygyrus' and organism[j] == 'haemonchus_contortus') or (organism[j] == 'heligmosomoides_polygyrus' and organism[i] == 'haemonchus_contortus'):

                if checkIfExist(2, organism[i], organism[j]):
                    j += 1
                    continue
                print(f'{organism[i]}, {organism[j]}')
                # set the results file path
                if settings.has_option('mode_2', 'file_1_path'):
                    settings.set('mode_2', 'file_1_path', output_mode1_path + '/results_' + organism[i] + '.csv')
                if settings.has_option('mode_2', 'file_2_path'):
                    settings.set('mode_2', 'file_2_path', output_mode1_path + '/results_' + organism[j] + '.csv')
                if settings.has_option('mode_2', 'organism_file_1'):
                    settings.set('mode_2', 'organism_file_1', organism[i])
                if settings.has_option('mode_2', 'organism_file_2'):
                    settings.set('mode_2', 'organism_file_2', organism[j])
                    with open('settings.ini', 'w') as settingFile:
                        settings.write(settingFile)
                # Call mode 2
                generate_mode2()
                j += 1
            i += 1

pipe()