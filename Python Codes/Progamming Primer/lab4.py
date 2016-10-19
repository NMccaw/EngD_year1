# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 13:31:36 2016

@author: nm2e16
"""
import math


def seq_sqrt(xs):
    """
    Returns the square root of each number of a list
    """
    num_list = []
    for xs_split in xs:
        print(xs)
        xs_num = int(xs_split)
        print(xs_num)
        xs_squrt = math.sqrt(xs_num)
        print(xs_squrt)
        num_list.append(xs_squrt)
    return num_list


def mean(xs):
    """
    Returns the average of a list of numbers
    """
    ave = 0
    for xs_split in xs:
        num = float(xs_split)
        print(xs_split)
        ave = ave+num
    average = ave/len(xs)
    return average


def wc(filename):
    """
    returns the word count of a text file
    """
    f = open(filename, 'rt')
    data = f.readlines()
    f.close()
    word_count_tot = 0
    for s in data:
        words = s.split()
        word_count = len(words)
        word_count_tot = word_count_tot+word_count
    return word_count_tot
