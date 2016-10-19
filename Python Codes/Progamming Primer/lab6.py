# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 12:06:57 2016

@author: nm2e16
"""


def eval_f(f, xs):
    """
    returns the results of a list after they have been throug a function
    """
    res_list = []
    for num in xs:
        #int_num = int(num)
        fun_num = f(num)
        res_list.append(fun_num)

    return res_list


def sum_f(f, xs):
    """
    returns the sum of a list after they have been through a function
    """
    sum_num = 0
    for num in xs:
        int_num = int(num)
        fun_num = f(int_num)
        sum_num = sum_num+fun_num
    return sum_num


def box_volume_UPS(a=13, b=11, c=2):
    """
    returns the volume of a box with standard measurements
    """
    vol = a*b*c
    return vol
