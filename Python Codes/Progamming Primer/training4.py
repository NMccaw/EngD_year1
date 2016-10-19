# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 11:34:24 2016

@author: nm2e16
"""
import urllib.request


def line_averages(filename):
    """
    Returns the average of each line
    """
    f = open(filename, 'rt')
    lines = f.readlines()
    f.close()
    my_list = []
    for line in lines:
        print(line)
        linesp = line.split(',')
        ave = 0
        for k in linesp:
            print(k)
            num = float(k)
            ave = num+ave
        average = ave/len(linesp)
        my_list.append(average)
        print(k)
    return my_list


def noaa_string():
    url = "http://tgftp.nws.noaa.gov/data/observations/metar/decoded/EGHI.TXT"
    noaa_data_string = urllib.request.urlopen(url).read()
    return noaa_data_string.decode("utf-8")


def noaa_temperature(s):
    """
    Returns the temperature from a report
    """
    print(s)
    lines = s.split('\n')
    for line in lines:
        print(line[0:3])
        if line[0:3] == 'Tem':
            imp_line = line
            print(imp_line)
            temp_str = imp_line.split(' ')
            temp_str2 = temp_str[3]
            numstr = temp_str2.strip('(')
            num2 = float(numstr)

    return num2
