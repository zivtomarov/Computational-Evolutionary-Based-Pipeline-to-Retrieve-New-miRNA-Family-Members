# coding: utf-8
import glob
import os
import configparser
import pandas as pd

organism = []

# read settings file
settings = configparser.ConfigParser()
settings._interpolation = configparser.ExtendedInterpolation()
settings.read('settings.ini')


# All Against All algorithm to gather all scores
def allAgainstAll():
    # go over all possible couples of organisms from mode2 results
    if settings.has_option('mode_4', 'dir_pairs_path'):
        dir_pairs_path = settings.get('mode_4', 'dir_pairs_path')

        if settings.has_option('mode_4', 'output_path'):
            output_path = settings.get('mode_4', 'output_path')

        if not os.path.exists(output_path):
            os.mkdir(output_path)

        # for each file inside directory create matrix
        filenames = glob.glob(dir_pairs_path + "/*.csv")


        list_all_pairs_scores = []


        for filename in filenames:
            print(filename)

            name_row = filename.split(dir_pairs_path + '/')[1].split(',')[0]
            name_column = filename.split(',')[1].split('.')[0]

            req_col = ['index_1', 'index_2', 'score_mature', 'score_hairpin'] #  'score_hairpin'
            df = pd.read_csv(filename, usecols=req_col)

            # run on all rows in pair_file
            for line in df.itertuples(index=False, name=None):
                # print(line)
                index_1 = line[0]
                index_2 = line[1]
                score_mature = line[2]
                hairpin_score = line[3]

                list_all_pairs_scores.append((name_row + "," + str(index_1), name_column + "," + str(index_2), score_mature, hairpin_score))


        if not os.path.exists(output_path + "/output"):
            os.mkdir(output_path + "/output")
        # write to file all pairs scores - mature and hairpin
        if not os.path.exists(output_path + "/output/pairs_score_mature.csv"):
            df_pairs_mature_score = pd.DataFrame(list_all_pairs_scores)
            df_pairs_mature_score.columns = ['organism1', 'organism2', 'mature_score', 'hairpin_score']
            df_pairs_mature_score.to_csv(output_path + "/output/pairs_scores.csv", index=False)



# allAgainstAll()
# print("Finished All-Against-All")

if settings.has_option('mode_4', 'output_path'):
    output_path = settings.get('mode_4', 'output_path')
df = pd.read_csv(output_path + "/output/pairs_score_mature.csv")
print("ye")
