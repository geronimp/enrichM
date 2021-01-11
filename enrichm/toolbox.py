#!/usr/bin/env python3
# pylint: disable=line-too-long
'''
Miscellaneous functions that dont necessarily belong anywhere but are still useful
'''
import logging
import subprocess
from itertools import islice

def list_splitter(input_list, chunk_number, chunk_max):
    """
    An iterator that separates a list into a number of smaller lists
    (chunk_number). Maximum size for the sub-lists can also be
    specified (chunk_max)
    
    :param input_list: List to be split into chunks (smaller lists)
    :type input_list: list
    :param chunk_number: Number of chunks to split list into
    :type chunk_number: int
    :param chunk_max: Maximum number of entries to have in each chunk.
    :type chunk_max: int
    :returns: A list, a chunk of input_list
    :rtype: list
    """
    list_size = float(len(input_list))
    chunk_size = int(round((list_size/chunk_number), 0))

    if chunk_size > chunk_max:
        chunk_size = chunk_max
    elif chunk_size < 1:
        chunk_size = list_size

    while list_size > 0:

        if len(input_list) <= chunk_size:
            yield input_list
            del input_list
        else:
            yield input_list[:chunk_size]
            del input_list[:chunk_size]

        try:
            list_size = len(input_list)
        except NameError:
            list_size = 0

def run_command(cmd):
    '''
    Runs a command using subprocess
    :param cmd: A string of a command to run.
    :type cmd: str
    :returns: None
    :rtype: bool
    '''
    logging.debug(cmd)
    subprocess.call(cmd, shell=True)
    logging.debug('Finished')

def get_present_annotations(input_dictionary):
    '''
    Takes an input dictionary and returns keys when they have a value of > 0
    :param input_dictionary: input dictionary, with arbitrary keys and int or float values
    :type input_dictionary: dict
    :returns: list of keys with value > 0
    :rtype: list
    '''
    output_list = list()

    for key, value in input_dictionary.items():

        if value > 0:
            output_list.append(key)

    return output_list


def window(seq, n=2):
    "Stolen directly from https://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator"
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def cluster(input_list, maxgap):
    '''Stolen from https://stackoverflow.com/questions/14783947/grouping-clustering-numbers-in-python'''
    input_list.sort()
    groups = [[input_list[0]]]
    for x in input_list[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups

def reverse_dictionary_of_lists(input_dictionary):
    output_dictionary = dict()

    for key, item_list in input_dictionary.items():
        item_list = list(item_list)

        for item in item_list:

            if item not in output_dictionary:
                output_dictionary[item] = list()

            output_dictionary[item].append(key)

    return output_dictionary