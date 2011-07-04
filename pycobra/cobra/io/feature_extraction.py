#cobra.io.feature_extraction.py
#Tools for parsing agilent feature extraction files.
from time import time
from rpy2.robjects import r
import rpy2.robjects.numpy2ri
from copy import deepcopy
import os, sys
from collections import defaultdict
from numpy import mean, array,std, log10
from cobra.tools import log_function
from cobra.stats.stats import combine_p_values, error_weighted
print "WARNING: cobra.io.feature_extraction is not ready for general use "
def parse_file(in_file, polarity=1, return_id='accession',
               normalization=None, lowess_parameter=0.33):
    """Extract data from feature extraction >=9.5 files.  Returns
    the average log ratios for each return_id.

    in_file: String.  Name of input file.

    polarity: 1 or -1.  Indicates whether to do red over green (1) or
    green over red.  If normalization isn't performed then the polarity
    is multiplied by the log ratio.

    return_id: 'accession' or 'probe'

    normalization: 'lowess' or None

    lowess_parameter: Float.  Smoothing parameter for lowess normalization

    """
    
    in_file_handle = open(in_file)
    the_header = in_file_handle.readline().rstrip('\r\n').split('\t')
    while not the_header[0] == 'FEATURES':
        the_header = in_file_handle.readline().rstrip('\r\n').split('\t')

    the_data = [x.rstrip('\r\n').split('\t')
                for x in in_file_handle.readlines() \
                if 'GE_BrightCorner' not in x and 'DarkCorner' not in x]
    in_file_handle.close()
    if return_id.lower() == 'accession':
        return_id = 'SystematicName'
    elif return_id.lower() == 'probe':
        return_id = 'ProbeName'
    else:
        return_id = 'SystematicName'
    gene_index = the_header.index(return_id)
    the_gene_dict = defaultdict(list)
    if not normalization:
        log_ratio_index = the_header.index('LogRatio')
        [the_gene_dict[x[gene_index]].append(float(x[log_ratio_index])*polarity)
         for x in the_data]
    elif normalization.lower() == 'lowess':
        red_index = the_header.index('rBGSubSignal')
        green_index = the_header.index('gBGSubSignal')
        gene_names = []
        red_signal = []
        green_signal = []
        [(gene_names.append(x[gene_index]),
          red_signal.append(float(x[red_index])),
          green_signal.append(float(x[green_index])))
         for x in the_data if float(x[red_index]) > 0 and float(x[green_index]) > 0]
        red_signal = array(red_signal)
        green_signal = array(green_signal)
        a = log10(red_signal*green_signal)/2. 
        if polarity == 1:
            m = log10(red_signal/green_signal)
        else:
            m = log10(green_signal/red_signal)
        lowess_fit = array(r.lowess(a, m, lowess_parameter)).T
        m = m - lowess_fit[:,1]
        #HERE: there's a problem with zipping
        [the_gene_dict[k].append(v) for k,v in zip(gene_names,list(m))]
    for k,v in the_gene_dict.items():
        the_gene_dict[k] = log10(mean(map(lambda x: 10**x, v)))
    return the_gene_dict

def parse_file_a(in_file, polarity = 1, quality_control=False):
    """

    quality_control:  Boolean.  Indicates whether to exclude probes that
    do not meet the quality control standards.
    
    TODO merge this with parse_file
    """
    with open(in_file) as in_file_handle:
        the_header = in_file_handle.readline().rstrip('\r\n').split('\t')
        while not the_header[0] == 'FEATURES':
            the_header = in_file_handle.readline().rstrip('\r\n').split('\t')

        the_data = [x.rstrip('\r\n').split('\t')
                     for x in in_file_handle.readlines()
                    if 'GE_BrightCorner' not in x and 'DarkCorner' not in x]

    probe_index = the_header.index('ProbeName')
    gene_index = the_header.index('SystematicName')
    column_to_index = {'log_ratio': the_header.index('LogRatio'),
                       'log_error': the_header.index('LogRatioError'),
                       'p_value': the_header.index('PValueLogRatio'),
                       'intensity_1': the_header.index('gProcessedSignal'),
                       'intensity_2': the_header.index('rProcessedSignal')}
    if quality_control:
        quality_control_indices = {}
        [quality_control_indices.update({the_header.index(k): 1})
         for k in map(lambda x: x + 'IsFound', ['r','g'])
         if k in the_header]
        [[quality_control_indices.update({the_header.index(k): 0})
         for k in map(lambda x: x + y, ['r','g'])
          if k in the_header]
         for y in ['IsSaturated', 'IsFeatNonUnifOL']]
        quality_control_indices[the_header.index('IsManualFlag')] = 0
        gene_dict = dict([(x[gene_index], defaultdict(list))
                          for x in the_data])
        for the_row in the_data:
            for k, v in quality_control_indices.items():
                if int(the_row[k]) != v:
                    continue
            the_dict = gene_dict[the_row[gene_index]]
            [the_dict[k].append(float(the_row[v]))
             for k, v in column_to_index.items()]
        [the_dict.pop(k)
         for k, v in the_dict.items()
         if len(v) == 0]
        
    else:
        gene_dict = dict([(x[gene_index], defaultdict(list))
                          for x in the_data])
        for the_row in the_data:
            the_dict = gene_dict[the_row[gene_index]]
            [the_dict[k].append(float(the_row[v]))
             for k, v in column_to_index.items()]

    [collapse_fields(v) for v in gene_dict.values()]
    return gene_dict



def parse_annotation(annotation_file, key_type='PrimaryAccession',
                     value_type='EntrezGeneID'):
    """Parse an agilent array annotation file into a dict.

    annotation_file: String.  Name of the annotation file.

    key_type: The column from the annotation_file to use as the key
    for the dict.

    value_type: The column from the annotation_file to use as the value
    for the dict.
    
    """
    annotation_dict = {}

    with open(annotation_file) as in_file:
        the_header = in_file.readline().rstrip('\r\n').split('\t')
        key_index = the_header.index(key_type)
        value_index = the_header.index(value_type)
        the_data = [x.rstrip('\r\n').split('\t') for x in in_file.readlines()]
        [annotation_dict.update({x[key_index]: int(x[value_index])})
         for x in the_data if len(x) >= value_index and x[value_index] != '']
    return annotation_dict


def collapse_fields(data_dict, quantitative_fields=['intensity_1',
                                                    'intensity_2'],
                    log_fields=['log_ratio', 'log_error'],
                    p_fields=['p_value'], log_base=10, error_weighting=True):
    """Collapses the rows for each field from a feature extraction element
    based on the data type.  Quantative fields are averaged, log_fields are
    converted to linear scale, averaged, and then log10 is taken. p-value
    fields are combined using the default method from cobra.stats.stats.
    combine_p_values.

    data_dict: A dictionary of numerical data values.

    quantitative_fields:  A list of the fields in data_dict that can be
    directly averaged.

    log_fields: A list of the fields in data_dict that are in log form
    and must be transformed before averaging.

    p_fields:  A list of the fields in data_dict that are p-values.

    log_base:  The base for the log transform.

    error_weighting:  Boolean.  If True return the error weighted-mean and error.
    Assumes the mean and std devs are the 1st two fields of log_fields.

    NOTE: This happens in place so the data_dict will be modified.

    """
    [data_dict.update({k: mean(data_dict[k])})
     for k in quantitative_fields]
    if error_weighting:
        log_fields = deepcopy(log_fields)
        mean_field = log_fields.pop(0)
        std_field = log_fields.pop(0)
        the_means = map(lambda x: log_base**x, data_dict[mean_field])
        the_stds = map(lambda x: log_base**x, data_dict[std_field])
        weighted_mean, weighted_std = error_weighted(the_means, the_stds)
        data_dict[mean_field] = log_function(weighted_mean, log_base)
        data_dict[std_field] = log_function(weighted_std, log_base)
    #TODO:  This needs to be updated for log_error
    [data_dict.update({k: log_function(mean(map(lambda x: log_base**x,
                                                data_dict[k])),
                                       log_base)})
     for k in log_fields]
    [data_dict.update({k: combine_p_values(data_dict[k])})
     for k in p_fields]

def combine_files(file_list, annotation_file, polarity_list=None,
                  print_time=False, quality_control=False):
    """Parse feature extraction files.  This function
    combines multiple technical replicates at the RNA level into a single
    experiment.  Typically multiple replicates will be used when dye-swapping
    is employed.

    file_list: A list of feature extraction file names.  Or a single file.
    A list is only provided if the values are to be combined across all files.

    annotation_file: The agilent annotation file corresponding to the feature
    extraction files in file_list.

    polarity_list:  A list of integers (1 or -1) indicating the polarity of
    the experiment.  1 indicates that the red channel is divided by the green
    channel and -1 indicates to divide the green channel by the red channel
    
    """
    if print_time:
        start_time = time()
    accession_to_entrez = parse_annotation(annotation_file)
    if print_time:
        print '%s %f'%('annotation time',
                       time() - start_time)
        start_time = time()

    parsed_files = []
    if not hasattr(file_list, '__iter__'):
        file_list = [file_list]
    if not hasattr(polarity_list, '__iter__'):
        if polarity_list is None:
            polarity_list = [1]*len(file_list)
        else:
            polarity_list = list(polarity_list)
    [parsed_files.append(parse_file_a(the_file,
                                     the_polarity,
                                      quality_control=quality_control))
     for the_file, the_polarity in zip(file_list,
                                       polarity_list)]
    if print_time:
        print '%s %f'%('parse time',
                       time() - start_time)
        start_time = time()

    #now combine measurements.  Assume that the files
    #all have the same ids in them.
    #Only look at items that are in the annotation file
    the_keys = set(parsed_files[0]).intersection(accession_to_entrez)
    combined_data = {}
    the_data_fields = parsed_files[0].values()[0].keys()
    #Merge the results for each field across files into a list
    for the_key in the_keys:
        tmp_dict = combined_data[the_key] = defaultdict(list)
        for data_dict in parsed_files:
            data_dict = data_dict[the_key]
            [tmp_dict[the_field].append(data_dict[the_field])
             for the_field in the_data_fields]
    if print_time:
        print '%s %f'%('combine time',
                       time() - start_time)
        start_time = time()

    #Collapse the lists for each field for each probe.
    [collapse_fields(v) for v in combined_data.values()]
    if print_time:
        print '%s %f'%('collapse time',
                       time() - start_time)
        start_time = time()
    
    #Now get the entrez ids for inserting into the database.
    the_entrez_ids = set([accession_to_entrez[the_key] for the_key in the_keys])
    #Now need to collapse to entrez gene ids
    entrez_data = dict([(k, defaultdict(list))
                        for k in the_entrez_ids])
    for the_id, the_dict in combined_data.items():
        tmp_dict = entrez_data[accession_to_entrez[the_id]]
        [tmp_dict[k].append(v)
         for k, v in the_dict.items()]
    #Now collapse the instances where the ids are have repeats.
    [collapse_fields(v) for v in entrez_data.values()]
    if print_time:
        print '%s %f'%('collapse to entrez and combine time',
                       time() - start_time)
    
    return entrez_data

