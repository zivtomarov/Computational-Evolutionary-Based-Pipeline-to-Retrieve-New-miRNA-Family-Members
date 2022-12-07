import os
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import pandas as pd
import configparser
import errors
import numpy as np
import numpy.matlib
from Bio import Align

# constant:
dirpath = os.getcwd()

# default parameters
file_1 = "file1.csv"
file_2 = "file2.csv"
min_alignment_score = 0
organism_file_1 = ""
organism_file_2 = ""

# [match, mismatch, start gap, gap extend]
# these are old params
# params_mature = [1, -2, -2, -2]
# params_hairpin = [2, -3, -5, -2]
params_mature = [1, -1, -1, -1]
params_hairpin = [1, -1, -1, -1]

# read settings file
settings = configparser.ConfigParser()
settings._interpolation = configparser.ExtendedInterpolation()
settings.read('settings.ini')

if not settings.has_section('mode_2'):
    print("mode_2 2 error: " + errors.errorMessage('section_not_exist'))

if settings.has_option('mode_2', 'file_1_path'):
    file_1 = settings.get('mode_2', 'file_1_path')
else:
    print("mode 2 error: " + errors.errorMessage('missing_input_file'))

if settings.has_option('mode_2', 'file_2_path'):
    file_2 = settings.get('mode_2', 'file_2_path')
else:
    print("mode 2 error: " + errors.errorMessage('missing_input_file'))

if settings.has_option('mode_2', 'organism_file_2'):
    organism_file_2 = settings.get('mode_2', 'organism_file_2')

if settings.has_option('mode_2', 'organism_file_1'):
    organism_file_1 = settings.get('mode_2', 'organism_file_1')

if settings.has_option('mode_2', 'output_path'):
    output_path = settings.get('mode_2', 'output_path')
if not os.path.isdir(output_path):
    os.mkdir(output_path)

if settings.has_option('mode_2', 'min_alignment_score'):
    min_alignment_score = float(settings.get('mode_2', 'min_alignment_score'))

if settings.has_option('mode_2', 'match_mature'):
    params_mature[0] = int(settings.get('mode_2', 'match_mature'))
if settings.has_option('mode_2', 'mismatch_mature'):
    params_mature[1] = int(settings.get('mode_2', 'mismatch_mature'))
if settings.has_option('mode_2', 'start_gap_mature'):
    params_mature[2] = int(settings.get('mode_2', 'start_gap_mature'))
if settings.has_option('mode_2', 'gap_extend_mature'):
    params_mature[3] = int(settings.get('mode_2', 'gap_extend_mature'))
if settings.has_option('mode_2', 'match_hairpin'):
    params_hairpin[0] = int(settings.get('mode_2', 'match_hairpin'))
if settings.has_option('mode_2', 'mismatch_hairpin'):
    params_hairpin[1] = int(settings.get('mode_2', 'mismatch_hairpin'))
if settings.has_option('mode_2', 'start_gap_hairpin'):
    params_hairpin[2] = int(settings.get('mode_2', 'start_gap_hairpin'))
if settings.has_option('mode_2', 'gap_extend_hairpin'):
    params_hairpin[3] = int(settings.get('mode_2', 'gap_extend_hairpin'))

if settings.has_option('mode_2', 'kth_score_mark'):
    kth_score_mark = int(settings.get('mode_2', 'kth_score_mark'))


if settings.has_option('mode_2', 'percentile_filter_percentage'):
    percentile_filter_percentage = float(settings.get('mode_2', 'percentile_filter_percentage'))

if settings.has_option('mode_2', 'value_filter'):
    value_filter = int(settings.get('mode_2', 'value_filter'))

if settings.has_option('mode_2', 'elbow_point'):
    elbow_point = settings.get('mode_2', 'elbow_point')

def get_max_score(alignments):
    max_score = 0
    max_alignment = ""
    for a in alignments:
        if max_score == 0:
            max_score = a[2]
            max_alignment = format_alignment(*a)
        else:
            if a[2] > max_score:
                max_score = a[2]
                max_alignment = format_alignment(*a)
    return [max_alignment, max_score]


def str_params(arr):
    res = "parameters:"
    index = 0
    for a in arr:
        if index == 0:
            res += str(a)
            index += 1
        else:
            res += "_"
            res += str(a)
    return res


def create_html_and_csv_file(res, res_file_name):
    # read settings file
    settings = configparser.ConfigParser()
    settings._interpolation = configparser.ExtendedInterpolation()
    settings.read('settings.ini')

    if settings.has_option('mode_2', 'organism_file_2'):
        organism_file_2 = settings.get('mode_2', 'organism_file_2')

    if settings.has_option('mode_2', 'organism_file_1'):
        organism_file_1 = settings.get('mode_2', 'organism_file_1')

    if settings.has_option('mode_2', 'match_mature'):
        match_mature = settings.get('mode_2', 'match_mature')
    if settings.has_option('mode_2', 'mismatch_mature'):
        mismatch_mature = settings.get('mode_2', 'mismatch_mature')
    if settings.has_option('mode_2', 'start_gap_mature'):
        start_gap_mature = settings.get('mode_2', 'start_gap_mature')
    if settings.has_option('mode_2', 'gap_extend_mature'):
        gap_extend_mature = settings.get('mode_2', 'gap_extend_mature')
    if settings.has_option('mode_2', 'match_hairpin'):
        match_hairpin = settings.get('mode_2', 'match_hairpin')
    if settings.has_option('mode_2', 'mismatch_hairpin'):
        mismatch_hairpin = settings.get('mode_2', 'mismatch_hairpin')
    if settings.has_option('mode_2', 'start_gap_hairpin'):
        start_gap_hairpin = settings.get('mode_2', 'start_gap_hairpin')
    if settings.has_option('mode_2', 'gap_extend_hairpin'):
        gap_extend_hairpin = settings.get('mode_2', 'gap_extend_hairpin')

    start_html = """<!DOCTYPE html>
    <html lang="en">
    <head>
        <style>
        table {
            text-align: center;
            font-size: 20px;
        }
        </style>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <meta http-equiv="X-UA-Compatible" content="ie=edge">
        <title>Document</title>
        <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.0.10/css/all.css" integrity="sha384-+d0P83n9kaQMCwj8F4RJB66tzIwOKmrdb46+porD/OvrJ+37WqIM7UoBtwHO6Nlg" crossorigin="anonymous">
        <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.0/css/bootstrap.min.css" integrity="sha384-9gVQ4dYFwwWSjIDZnLEWnxCjeSWFphJiwGPXr1jddIhOegiu1FwO5qRGvFXOdJZ4" crossorigin="anonymous">    
    </head>
    <body>
        <div class="container" style="max-width:1500px;">
            <br>
            <h1><i class="fas fa-dna"></i> &nbsp; Alignment Results</h1>
            <br>
            <br>
            <h4> &nbsp; match mature: """ + match_mature + """, mismatch_mature: """ + mismatch_mature + """, start_gap_mature: """ + start_gap_mature + """, gap_extend_mature: """ + gap_extend_mature + """</h4>
            <br>
            <table class="table table-striped">
                <thead>
                 <tr>
                        <th></th>
                        <th colspan=3>""" + organism_file_1 + """</th>
                        <th colspan=3>""" + organism_file_2 + """</th>
                        <th colspan =3>mature</th>
                    </tr>
                    <tr>
                        <th></th>
                        <th>index</th>
                        <th>type</th>
                        <th>found in DB</th>
                        <th>index</th>
                        <th>type</th>
                        <th>found in DB</th>
                        <th>alignment</th>
                        <th>score mature</th>
                        <th>score hairpin</th>

                    </tr>
                </thead>
                <tbody>"""

    end_html = """</tbody>
            </table>
        </div>
    </body>
    </html>"""


    start_csv = "mature_1,hairpin_1,index_1,type_1,found_in_DB_1,mature_2,hairpin_2,index_2,type_2,found_in_DB_2,score_mature,score_hairpin\n"

    try:
        write_file_html = open(output_path + '/' + res_file_name + ".html", "w")
        write_file_csv = open(output_path + '/' + res_file_name + ".csv", "w")
        write_file_html.write(start_html)
        write_file_csv.write(start_csv)
        i = 1
        for row in res:
            allignment_mature = row[10].split("S")[0].rstrip()
            # allignment_mature = row[10]
            #
            # # allignment_hairpin = row[13].split("S")[0].rstrip() # removed from mode2
            #
            # # allignment_hairpin = row[13]
            write_file_html.write("<tr>")
            write_file_html.write("<td>" + str(i) + "</td>\n")
            write_file_html.write("<td>" + str(row[2]) + "</td>\n")
            write_file_html.write("<td>" + str(row[3]) + "</td>\n")
            write_file_html.write("<td>" + str(row[4]).replace('nan', '') + "</td>\n")
            write_file_html.write("<td>" + str(row[7]) + "</td>\n")
            write_file_html.write("<td>" + str(row[8]) + "</td>\n")
            write_file_html.write("<td>" + str(row[9]).replace('nan', '') + "</td>\n")
            write_file_html.write("<td><pre>" + allignment_mature + "</pre></td>\n")
            write_file_html.write("<td>" + str(row[11]) + "</td>\n")
            write_file_html.write("<td>" + str(row[12]) + "</td>\n")
            # write_file_html.write("<td>" + str(row[12]) + "</td>\n")
            # write_file_html.write("<td><pre>" + allignment_hairpin + "</pre></td>\n")
            # write_file_html.write("<td>" + str(row[14]) + "</td>\n")
            # write_file_html.write("<td>" + str(row[15]) + "</td>\n")
            write_file_html.write("</tr>")


            # write to csv all beside the alignments themselves
            join_by = ","
            # csv = join_by.join(row).replace(','+str(row[10]),'').replace(','+str(row[13]),'')

            # csv = join_by.join( # removed from mode2
            #     [str(row[0]), str(row[1]), str(row[2]), str(row[3]), str(row[4]).replace('nan', ''), str(row[5]),
            #      str(row[6]), str(row[7]), str(row[8]), str(row[9]).replace('nan', ''), str(row[11]), str(row[12]),
            #      str(row[14]), str(row[15])]) + "\n"
            csv = join_by.join(
                [str(row[0]), str(row[1]), str(row[2]), str(row[3]), str(row[4]).replace('nan', ''), str(row[5]),
                 str(row[6]), str(row[7]), str(row[8]), str(row[9]).replace('nan', ''), str(row[11]), str(row[12])]) + "\n"
            write_file_csv.write(csv)

            i = i + 1

        write_file_html.write(end_html)
        write_file_html.close()
        write_file_csv.close()
    finally:
        write_file_html.close()
        write_file_csv.close()


def find_cut_off(array, sum_array, percentage):
    cut_off = sum_array * percentage

    # print("sum_array : ", sum_array)
    # print("percentage : ", percentage)
    # print("cut_off : ", cut_off)

    sum_num = 0
    index = len(array) - 1
    while index > 0:
        if array[index] + sum_num > cut_off:
            break
        else:
            sum_num += array[index]
            index = index - 1
    return index # + 1


# gets 2 files generated from pipeline and compare all against all
# parameters are different for hairpin/mature
# mode can be hairpin compare or mature compare
def align_files(file_1, file_2, name_organism_1, name_organism_2, percentile_filter_percentage, value_filter, elbow_point):
    # read settings file
    settings = configparser.ConfigParser()
    settings._interpolation = configparser.ExtendedInterpolation()
    settings.read('settings.ini')

    if settings.has_option('mode_2', 'threshold'):
        threshold = int(settings.get('mode_2', 'threshold'))

    if settings.has_option('mode_2', 'hairpin_filter'):
        hairpin_filter = int(settings.get('mode_2', 'hairpin_filter'))

    if settings.has_option('mode_2', 'kth_score_mark'):
        kth_score_mark = int(settings.get('mode_2', 'kth_score_mark'))

    if settings.has_option('mode_2', 'match_mature'):
        match_mature = int(settings.get('mode_2', 'match_mature'))
    if settings.has_option('mode_2', 'mismatch_mature'):
        mismatch_mature = int(settings.get('mode_2', 'mismatch_mature'))
    if settings.has_option('mode_2', 'gap_score_mature'):
        gap_score_mature = int(settings.get('mode_2', 'gap_score_mature'))

    temp_results = []

    join_by = ","
    res_file_name = join_by.join([name_organism_1, name_organism_2])

    # create array for each score and increase by 1 (counter) when score is found
    array = [0] * 23  # for scores 0 -> 22


    # req_col = ['Unnamed: 0', 'mature_3p', 'mature_5p', 'seq_cut_tails', 'found_in_DB_id']

    req_col = ['Unnamed: 0', 'Mature', '3p/5p', 'Hairpin_seq_trimmed', 'Found_in_DB_id']

    # loop over all seq in files
    df1 = pd.read_csv(file_1, usecols=req_col)  # df1 = pd.read_csv(dirpath + "/" + file_1)
    df2 = pd.read_csv(file_2, usecols=req_col)
    i = 0
    print("df1 len : ", len(df1))
    print("df2 len : ", len(df2))

    # with open('temp_csv_file.csv', "w") as f:
    #     writer = csv.writer(f)
    #     writer.writerow(['index_1', 'index_2', 'alignment_mature'])

    # configure the aligner object
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'

    # aligner.match = 1.0
    # aligner.mismatch = -1.0
    # aligner.gap_score = -1.0

    aligner.match = match_mature
    aligner.mismatch = mismatch_mature
    aligner.gap_score = gap_score_mature

    for line_f1 in df1.itertuples(name=None):
        if i % 100 == 0:
            print(i)
        i += 1
        for line_f2 in df2.itertuples(name=None):

            # at the pairs of the same organism don't take organism x with organism x
            if name_organism_1 == name_organism_2 and line_f1[1] == line_f2[1]:
                continue
            if not pd.isnull(line_f1[2]):
                mature_1 = line_f1[2]
            else:
                mature_1 = line_f1[3]
            if not pd.isnull(line_f2[2]):
                mature_2 = line_f2[2]
            else:
                mature_2 = line_f2[3]

            index_1 = line_f1[1]  # index_id
            index_2 = line_f2[1]  # index_id


            # mature alignment
            # alignment_mature_pairwise = pairwise2.align.globalms(mature_1, mature_2, params_mature[0], params_mature[1],
            #                                              params_mature[2], params_mature[3], one_alignment_only=True)
            #
            # alignment_mature = format_alignment(*alignment_mature_pairwise[0])
            # score_mature = int(alignment_mature_pairwise[0][2])


            alignment_mature = aligner.align(mature_1, mature_2)
            score_mature = alignment_mature[0].score


            if score_mature >= threshold:

                newline = [index_1, index_2, str(alignment_mature[0]), score_mature]
                temp_results.append(newline)

                # writer.writerow(newline)

            if score_mature <= 0:
                array[0] += 1

            else:
                array[int(score_mature)] += 1

        # drop the candidates already paired so won't be a->b and b->a for same organism
        if name_organism_1 == name_organism_2:
            df2.drop(df2.index[0], inplace=True)


    print("finished long loop")

    # algorithm for percentile with percentile only

    if kth_score_mark > 0:
        ind = [index for index, item in enumerate(array) if item != 0][-3]
        print("cutoff score {} for k largest number".format(ind))
        results = []
        for i in temp_results:
            if i[3] >= ind:  # results by mature score
                results.append(i)
    else:

        if value_filter == 0:
            sum_array = sum(array)
            index_percentile = find_cut_off(array, sum_array, percentile_filter_percentage)  # final cut off value represents the 95% top scores (95th percentile)
            print("cutoff score {} for percentile {}".format(index_percentile, percentile_filter_percentage))


            results = []
            for i in temp_results:
                if i[3] >= index_percentile:  # results by mature score
                    results.append(i)
        else:  # when user wants threshold of specific value
            results = []
            for i in temp_results:
                if i[3] >= value_filter:  # results by mature score
                    results.append(i)

    # another filtering using elbow point
    if elbow_point == 'True' and index_percentile < len(array) - 1:
        # print("potential number of candidates : ", len(temp_results))
        # create_html_and_csv_file(temp_results, res_file_name)

        # algorithm for elbow point
        curve = []
        for sublist in results:
            curve.append(sublist[3])
        curve.sort(reverse=True)
        nPoints = len(curve)
        allCoord = np.vstack((range(nPoints), curve)).T
        np.array([range(nPoints), curve])
        # pull out first point
        firstPoint = allCoord[0]
        # get vector between first and last point - this is the line
        lineVec = allCoord[-1] - allCoord[0]
        # normalize the line vector
        lineVecNorm = lineVec / np.sqrt(np.sum(lineVec ** 2))
        # find the distance from each point to the line:
        # vector between all points and first point
        vecFromFirst = allCoord - firstPoint
        # To calculate the distance to the line, we split vecFromFirst into two
        # components, one that is parallel to the line and one that is perpendicular
        # Then, we take the norm of the part that is perpendicular to the line and get the distance.
        # We find the vector parallel to the line by projecting vecFromFirst onto
        # the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
        # We project vecFromFirst by taking the scalar product of the vector with
        # the unit vector that points in the direction of the line (this gives us
        # the length of the projection of vecFromFirst onto the line). If we
        # multiply the scalar product by the unit vector, we have vecFromFirstParallel
        scalarProduct = np.sum(vecFromFirst * np.matlib.repmat(lineVecNorm, nPoints, 1), axis=1)
        vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
        vecToLine = vecFromFirst - vecFromFirstParallel
        # distance to line is the norm of vecToLine
        distToLine = np.sqrt(np.sum(vecToLine ** 2, axis=1))
        # find the maximum
        idxOfBestPoint = np.argmax(distToLine)

        score = curve[idxOfBestPoint-1]
        elbow_point_results = []
        for i in results:
            if i[3] >= score:  # results by mature score
                elbow_point_results.append(i)

        print("potential number of candidates : ", sum(array))
        print("number of candidates filtered by percentile : ", len(results))
        print("number of candidates filtered by elbow: ", len(elbow_point_results))
        print("cutoff score {} for elbow point".format(score))


        # gather all data from dataframes in the memory - df1 & df2
        final_results = []

        for i in elbow_point_results:
            organism1_info = df1[df1['Unnamed: 0'].eq(i[0])].values[0]
            mature_1 = organism1_info[1]
            type_1 = organism1_info[2]

            # if not pd.isnull(organism1_info[2]):
            #     mature_1 = organism1_info[2]  # mature_3p
            #     type_1 = '3p'
            # else:
            #     mature_1 = organism1_info[4]  # mature_5p
            #     type_1 = '5p'
            hairpin_1 = organism1_info[3]
            index_1 = organism1_info[0]
            found_in_DB_1 = organism1_info[4]

            if name_organism_1 == name_organism_2:
                organism2_info = df1[df1['Unnamed: 0'].eq(i[1])].values[0]
            else:
                organism2_info = df2[df2['Unnamed: 0'].eq(i[1])].values[0]

            mature_2 = organism2_info[1]
            type_2 = organism2_info[2]

            # if not pd.isnull(organism2_info[2]):
            #     mature_2 = organism2_info[2]  # mature_3p
            #     type_2 = '3p'
            # else:
            #     mature_2 = organism2_info[4]  # mature_5p
            #     type_2 = '5p'
            hairpin_2 = organism2_info[3]
            index_2 = organism2_info[0]
            found_in_DB_2 = organism2_info[4]

            alignment_mature = i[2]
            score_mature = i[3]

            alignment_hairpin = aligner.align(hairpin_1, hairpin_2)
            score_hairpin = alignment_hairpin[0].score

            newline = [mature_1, hairpin_1, index_1, type_1, found_in_DB_1, mature_2, hairpin_2, index_2, type_2,
                       found_in_DB_2, alignment_mature, score_mature, score_hairpin]
            final_results.append(newline)


        create_html_and_csv_file(final_results, res_file_name)

    else:
        # results = sorted(results, key=itemgetter(11), reverse=True)  # sort results by homology
        print("potential number of candidates : ", sum(array))
        print("number of candidates filtered: ", len(results))

        # gather all data from dataframes in the memory - df1 & df2
        final_results = []
        hairpin_scores = []
        for i in results:
            organism1_info = df1[df1['Unnamed: 0'].eq(i[0])].values[0]
            # organism1_info = df1.iloc[i[0]].values
            mature_1 = organism1_info[1]
            type_1 = organism1_info[2]
            # if not pd.isnull(organism1_info[2]):
            #     mature_1 = organism1_info[2]  # mature_3p
            #     type_1 = '3p'
            # else:
            #     mature_1 = organism1_info[1]  # mature_5p
            #     type_1 = '5p'
            hairpin_1 = organism1_info[3]
            index_1 = organism1_info[0]
            found_in_DB_1 = organism1_info[4]

            if name_organism_1 == name_organism_2:
                organism2_info = df1[df1['Unnamed: 0'].eq(i[1])].values[0]
            else:
                organism2_info = df2[df2['Unnamed: 0'].eq(i[1])].values[0]

            mature_2 = organism2_info[1]
            type_2 = organism2_info[2]

            # if not pd.isnull(organism2_info[2]):
            #     mature_2 = organism2_info[2]  # mature_3p
            #     type_2 = '3p'
            # else:
            #     mature_2 = organism2_info[1]  # mature_5p
            #     type_2 = '5p'
            hairpin_2 = organism2_info[3]
            index_2 = organism2_info[0]
            found_in_DB_2 = organism2_info[4]

            alignment_mature = i[2]
            score_mature = i[3]

            alignment_hairpin = aligner.align(hairpin_1, hairpin_2)
            score_hairpin = alignment_hairpin[0].score

            if hairpin_filter > 0:
                alignment_hairpin = aligner.align(hairpin_1, hairpin_2)
                score_hairpin = alignment_hairpin[0].score
                hairpin_scores.append(score_hairpin)

                if score_hairpin >= hairpin_filter:
                    newline = [mature_1, hairpin_1, index_1, type_1, found_in_DB_1, mature_2, hairpin_2, index_2,
                               type_2,
                               found_in_DB_2, alignment_mature, score_mature]
                    final_results.append(newline)


            else:
                newline = [mature_1, hairpin_1, index_1, type_1, found_in_DB_1, mature_2, hairpin_2, index_2, type_2,
                           found_in_DB_2, alignment_mature, score_mature, score_hairpin]
                final_results.append(newline)


        print("number of candidates after hairpin filtered: ", len(final_results))
        create_html_and_csv_file(final_results, res_file_name)





# start = time.clock()
# align_files(file_1, file_2, organism_file_1, organism_file_2, percentile_filter_percentage, value_filter, elbow_point)
# elapsed = (time.clock() - start)
# print("Program executed in " + str(elapsed))
# print("finish running mode 1")



def find_number_of_orgs(file_1, file_2, organism_file_1, organism_file_2):

    req_col = ['index_1', 'index_2']
    df_mode_2 = pd.read_csv(output_path + '/' + organism_file_1 + ',' + organism_file_2 + '.csv', usecols=req_col)


    req_col = ['Unnamed: 0', 'Found_in_DB_id']
    df = pd.read_csv(file_1, usecols=req_col)  # df1 = pd.read_csv(dirpath + "/" + file_1)
    counter_org_1 = 0
    counter_pos_org_1 = 0
    print("-----------------------------------------------")
    print("running over :", organism_file_1)
    print("candidates :", len(df))
    unique_set = df_mode_2['index_1'].unique()
    for line in df.itertuples(name=None):
        if line[1] in unique_set:
            counter_org_1 += 1
            if not pd.isnull(line[2]):
                counter_pos_org_1 += 1


    df = pd.read_csv(file_2, usecols=req_col)  # df1 = pd.read_csv(dirpath + "/" + file_1)

    print("running over :", organism_file_2)
    print("candidates :", len(df))
    counter_org_2 = 0
    counter_pos_org_2 = 0
    unique_set = df_mode_2['index_2'].unique()
    for line in df.itertuples(name=None):
        if line[1] in unique_set:
            counter_org_2 += 1
            if not pd.isnull(line[2]):
                counter_pos_org_2 += 1

    print("organism : {}, numer of candidates : {}, number of positive {}".format(organism_file_1,counter_org_1,counter_pos_org_1))
    print("organism : {}, numer of candidates : {}, number of positive {}".format(organism_file_2,counter_org_2,counter_pos_org_2))


# find_number_of_orgs(file_1, file_2, organism_file_1, organism_file_2)


# Count 5p and 3p from csv output files from mode 1
# read settings file
settings = configparser.ConfigParser()
settings._interpolation = configparser.ExtendedInterpolation()
settings.read('settings.ini')

def counter_report(csv_file):
    df = pd.read_csv(csv_file)
    ser = df['3p/5p'].value_counts()
    return {'3p': ser['3p'], '5p': ser['5p']}


if settings.has_option('mode_4', 'dir_pairs_path'):
    dir_result_mode2_path = settings.get('mode_4', 'dir_pairs_path')

# go over all files in mode1 to get the 3p and 5p candidates
for dirpath, dirs, files in os.walk(dir_result_mode2_path):
    for fname in files:
        file_name = os.path.join(dirpath,fname)
        if file_name.endswith('.csv'):
            # file_name = 'result_files_mode1_mir35/results_caenorhabditis_briggsae.csv'
            print('file_name {} report counter {}'.format(file_name, counter_report(file_name)))
            # print("*********************************")
