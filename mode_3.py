# coding: utf-8
import glob
import os
import configparser
import pandas as pd
import pybedtools as pybedtools


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





allAgainstAll()

print("Finished All-Against-All")

# if settings.has_option('mode_4', 'output_path'):
#     output_path = settings.get('mode_4', 'output_path')
# df = pd.read_csv(output_path + "/output/pairs_score_mature.csv")
# print("ye")



# def find_in_gff(gff_file, candidate_file, organism):
#
#     df_candidates = pd.read_csv(candidate_file)
#     # get only the relevant organism
#     df_candidates = df_candidates[df_candidates['Organism'].str.contains(organism)]
#
#     # data = {'chromosome': df_candidates['Chromosome'], 'source': "", 'type': "", 'start': df_candidates['Start_hairpin'],
#     #         'end': df_candidates['End_hairpin'], 'score': 0, 'strand': df_candidates['Strand'], 'phase': "", 'index': df_candidates['Organism']}
#     data = {'chromosome': df_candidates['Chromosome'], 'start': df_candidates['Start_hairpin'], 'end': df_candidates['End_hairpin'],
#             'strand': df_candidates['Strand'], 'index': df_candidates['Organism']}
#     df_candidates = pd.DataFrame(data, columns=['chromosome', 'start', 'end', 'strand', 'index'])
#     candidates = pybedtools.BedTool.from_dataframe(df_candidates)
#
#     candidates.saveas('final_candidates.bed')
#
#     print('## Candidates ##')
#     print(candidates)
#
#     # gff_bed = pybedtools.BedTool(gff_file)
#
#     # Parsing the GFF file to dataframe to bed file
#     col_list = ["chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
#     df_gff_init = pd.read_csv(gff_file, nrows=1000, compression='gzip', sep="\t", header=None, names=col_list, index_col=False)
#     df_gff_init = df_gff_init[~df_gff_init['chr'].str.contains('#')].dropna()
#
#     data = {'chromosome': df_gff_init['chr'], 'start': df_gff_init['start'], 'end': df_gff_init['end'],
#             'strand': df_gff_init['strand'], 'attributes': df_gff_init['attributes']}
#     df_gff = pd.DataFrame(data, columns=['chromosome', 'start', 'end', 'strand', 'attributes'])
#     df_gff['start'] = df_gff['start'].astype(int)
#     df_gff['end'] = df_gff['end'].astype(int)
#
#     organism_gff = pybedtools.BedTool.from_dataframe(df_gff)
#     # organism_gff = organism_gff.sort()
#     organism_gff.saveas('elegans_gff.bed')
#
#     print('## GFF ##')
#     print(organism_gff)
#
#     print(f"candidates - {candidates.fn}, organism_gff - {organism_gff.fn}")
#
#     # get the closest feature in 'other.bed' on the same strand
#     result = candidates.sort().closest(organism_gff.sort(), s=True, d=True) # D='b'
#     # results = dict([ (r.name, int(r[len(r.fields)-1])) for r in result])
#
#     print('## Results ##')
#     print("candidates in gff: ", len(result))
#     print(result)
#
#     result.saveas('candidates_in_gff.bed')


# gff_file = 'caenorhabditis_elegans.PRJNA13758.WBPS17.annotations.gff3.gz'
# organism = 'caenorhabditis_elegans'
# candidate_file = 'Final_candidates_mir35.csv'
# # candidates_file = 'intersection_full_overlap.bed'
# find_in_gff(gff_file, candidate_file, organism)



# check = pybedtools.BedTool('c_elegans_WS199_shortened_gff.txt')
#
# col_list = ["chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
# df_gff_init = pd.read_csv('c_elegans_WS199_shortened_gff.txt', sep="\t", header=None, names=col_list, index_col=False)
# # data = {'chromosome': df_gff_init['chr'], 'start': df_gff_init['start'], 'end': df_gff_init['end'],
# #         'name': "", 'score': 0, 'strand': df_gff_init['strand'], 'bio_type': df_gff_init['type']}
# # df_gff = pd.DataFrame(data, columns=['chromosome', 'start', 'end', 'name', 'score', 'strand', 'bio_type'])
#
# data = {'chromosome': df_gff_init['chr'], 'start': df_gff_init['start'], 'end': df_gff_init['end'],
#         'strand': df_gff_init['strand'], 'attributes': df_gff_init['attributes']}
# df_gff = pd.DataFrame(data, columns=['chromosome', 'start', 'end', 'strand', 'attributes'])
#
# check_pybedtools = pybedtools.BedTool.from_dataframe(df_gff)
# check_pybedtools.saveas('check_pybedtools.bed')

#



