# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 16:35:50 2016

@author: nm2e16
"""
import sys


def count_sub_in_file(filename, s):
    """
    return the word count of a file
    """

    f = open(filename, 'rt')

    data = f.readlines()
    f.close

    count = 0
    for words in data:
        print(words)

        word_split = words.split()
        print(word_split)
        for d in word_split:
            if d == s:
                count = count+1

    return count


def count_sub_in_file(filename, s):
    """
    return the word count of a file
    """
    try:
        f = open(filename, 'rt')
    except:
        return -1
    data = f.readlines()
    f.close

    count = 0
    for words in data:
        print(words)

        word_split = words.split()
        print(word_split)
        for d in word_split:
            if d == s:
                count = count+1

    return count


def count_vowels(s):
    """
    returns the vowel count
    """
    words = s.split()
    print(words)
    count = 0
    for d in words:

        for b in d:
            print(b)
            if b in ['a', 'e', 'i', 'o', 'u', 'A', 'E', 'I', 'O', 'U']:
                count = count+1
    return count
