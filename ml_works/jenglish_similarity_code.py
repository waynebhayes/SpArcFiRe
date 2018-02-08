import csv
import numpy as np
from numpy import isnan
from numpy import isinf
from scipy.stats.mstats import zscore
from scipy import spatial


def name_list_generator(file_name):
    datafile = open(file_name, 'r')
    datareader = csv.reader(datafile)
    data_1 = []
    for row in datareader:
        data_1.append(row)

    return data_1


def import_smoothed_file():
    import h5py
    h5f = h5py.File(r"dataV2.jld", 'r')
    dataFilteredSTD = np.transpose(h5f['dataFilteredSTD'][:])
    dataFilteredNames = np.transpose(h5f['dataFilteredNames'][:])
    bigYArray = np.transpose(h5f['bigYArray'][:])
    dataNames = np.transpose(h5f['namesArray'][:])
    bigYNames = np.transpose(h5f['bigYNames'][:])
    return dataFilteredSTD, dataFilteredNames, bigYArray, bigYNames, dataNames


def generate_remove_list(values):
    remove_list = []
    for x, testline in enumerate(values):
        try:
            if hasattr(np.float32(float(testline)), '__iter__'):
                remove_list.append(x)
        except:
            remove_list.append(x)

    return remove_list


def remove_strings(values, remove_list):
    return [val for x, val in enumerate(values) if x not in remove_list]


def treat_strings(values, remove_list, pa_list):
    removed = remove_strings(values, remove_list)
    return [float(val) if x in pa_list else val for x, val in enumerate(removed)]


def treat_strings_2d(values, remove_list, pa_list):
    data = []
    for val in values:
        data.append(treat_strings(val, remove_list, pa_list))
    return data


def smooth_test_features(test_feature_data):
    data = []
    for i, test_line in enumerate(test_feature_data):
        dat = []
        for test in test_line:
            if test.lower() == "true":
                dat.append(1)
            elif test.lower() == "false":
                dat.append(0)
            else:
                ft = np.float32(test)
                if isnan(ft):
                    dat.append(0)
                elif isinf(ft):
                    dat.append(3.40282e+38)
                else:
                    dat.append(ft)
        data.append(dat)
    return data


def standardize_values(values):
    vals = np.asarray(values, dtype=np.float32)
    for val in range(len(vals[0])):
        vals[:, val] = zscore(vals[:, val])
    return list(vals)


def data_maker(illustris_file_name, row_length=156, label_index=-1):
    if ".csv" in illustris_file_name:
        ill_names, ill_values, ill_labels, ill_names_0 = ill_camera_angle_data_generator_csv(illustris_file_name,
                                                                                             row_length, label_index)
    else:
        ill_names, ill_values, ill_labels, ill_names_0 = ill_camera_angle_data_generator_tsv(illustris_file_name,
                                                                                             row_length, label_index)
    remove_list = generate_remove_list(ill_values[0])
    treated_ill_labels = remove_strings(ill_labels, remove_list)
    pa_list = [x for x, label in enumerate(treated_ill_labels) if "pa_" in treated_ill_labels]
    abs_ill_values = treat_strings_2d(ill_values, remove_list, pa_list)
    abs_ill_values = smooth_test_features(abs_ill_values)
    # treated_ill_values = standardize_values(abs_ill_values)
    # return ill_names, treated_ill_values, treated_ill_labels, ill_names_0

    return ill_names, abs_ill_values, treated_ill_labels, ill_names_0


def strip_2d(two_d_list):
    return [list(map(str.strip, line)) for line in two_d_list]


def ill_camera_angle_data_generator_tsv(file_name, row_length=156, label_index=-1):
    datafile = open(file_name, 'r')
    datareader = csv.reader(datafile, delimiter="\t")
    data_1 = []
    for i, row in enumerate(datareader):
        if i == 0:
            data_1.append([r for r in row])
        elif float(row[143]) + float(row[144]) > 0.95:
            data_1.append([r for r in row])
    # print(len(data_1))
    # print(len(data_1[int(len(data_1)/2)]))
    data = [dat for dat in data_1 if len(dat) == row_length]
    data_0 = [dat for dat in data_1 if len(dat) != row_length]
    names = [dat[0] for dat in data]
    names_0 = [dat[0] for dat in data_0]
    if label_index == -1:
        values = [dat for dat in data[:-1]]
    else:
        values = [dat for dat in data[1:]]

    labels = data[label_index]
    return list(map(str.strip, names)), strip_2d(values), list(map(str.strip, labels)), list(map(str.strip, names_0))


def ill_camera_angle_data_generator_csv(file_name, row_length=156, label_index=-1):
    datafile = open(file_name, 'r')
    datareader = csv.reader(datafile)
    data_1 = []
    for row in datareader:
        data_1.append([r for r in row])
    # print(len(data_1[0]))
        data = [dat for dat in data_1 if len(dat) == row_length]
    data_0 = [dat for dat in data_1 if len(dat) != row_length]
    names = [dat[0] for dat in data[:-1]]
    names_0 = [dat[0] for dat in data_0]
    values = [dat for dat in data[:-1]]

    labels = data[label_index]
    return list(map(str.strip, names)), strip_2d(values), list(map(str.strip, labels)), list(map(str.strip, names_0))


def cosine_finder(gz, ill):
    return 1 - spatial.distance.cosine(gz, ill)


def min_max_maker(gz_co_2, ill_co_2):
    data = []
    for ill in ill_co_2:
        min = (0,0,999999)
        max = (0,0,-9999999)
        for gz in gz_co_2:
            if(gz[0] != ill[0]):
                gz_np = np.array(gz[1:], dtype=float)
                ill_np = np.array(ill[1:], dtype=float)
                point = (gz[0], ill[0], cosine_finder(gz_np, ill_np))
                if point[2] < min[2]:
                    min = point
                if point[2] > max[2]:
                    max = point
        data.append([min, max])
    return data


def distance_dict_maker(names, values):
    cosine_dict = {}
    for i in range(len(values)):
        cosine_values = []
        for y in range(len(values)):
            if i != y:
                cos = cosine_finder(values[i], values[y])
                if cos > 0.99:
                    cosine_values.append((names[y], cos))
        cosine_dict[names[i]] = sorted(cosine_values, key=lambda x: x[1], reverse=True)
        print(names[i], len(cosine_values), cosine_values[:4])
    return cosine_dict


def index_array_from_intersect_and_labels(intersect, labels):
    return [labels.index(intersect[i]) for i in range(len(intersect))]


def intersect_list_maker(intersect_list):
    intersect_output = intersect_list[0]
    for x in range(len(intersect_list) - 1):
        intersect_output = np.intersect1d(intersect_output, intersect_list[x+1])
    return intersect_output


def feature_data_generator(data, index_array):
    train_feature_data = []
    for d in data:
        train_feature_data.append([d[i] for i in index_array])
    return train_feature_data


def illustris_data_maker(illustris_file_name):
    if ".csv" in illustris_file_name:
        ill_names, ill_values, ill_labels, ill_names_0 = ill_camera_angle_data_generator_csv(illustris_file_name)
    else:
        ill_names, ill_values, ill_labels, ill_names_0 = ill_camera_angle_data_generator_tsv(illustris_file_name)
    remove_list = generate_remove_list(ill_values[0])
    treated_ill_labels = remove_strings(ill_labels, remove_list)
    pa_list = [x for x, label in enumerate(treated_ill_labels) if "pa_" in treated_ill_labels]
    abs_ill_values = treat_strings_2d(ill_values, remove_list, pa_list)
    abs_ill_values = smooth_test_features(abs_ill_values)
    treated_ill_values = standardize_values(abs_ill_values)
    return ill_names, treated_ill_values, treated_ill_labels, ill_names_0


def main():
    dataFilteredSTD, dataFilteredNames, bigYArray, bigYNames, dataNames = import_smoothed_file()
    ill_data_name = "ILLUSTRIS-SpArcFiRe-all/g/galaxy.csv"
    ill_names, ill_values, ill_labels, ill_names_0 = illustris_data_maker(ill_data_name)
    print("ill_done")
    gz_data_name = "allColumn-complete-MR-22.25-vs-z.085-all-pr90.gt.6.tsv"
    gz_names, gz_values, gz_labels, gz_names_0 = data_maker(gz_data_name, row_length=195, label_index=0)
    print("gz_done")
    intersects = intersect_list_maker([ill_labels, gz_labels, dataFilteredNames])
    print("intersects done:", str(len(intersects)))
    gz_indexes = index_array_from_intersect_and_labels(intersects, list(gz_labels))
    print("gz index done:", str(len(gz_indexes)))
    gz_feature_data = feature_data_generator(gz_values, gz_indexes)
    print("gz features data done:", str(len(gz_feature_data)))
    print([gz_labels[i] for i in gz_indexes])

    distance_dict = distance_dict_maker(gz_names[1:], gz_feature_data)
    print(distance_dict[gz_names[1]])


if __name__ == '__main__':
    main()
