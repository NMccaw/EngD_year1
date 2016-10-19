# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 11:44:18 2016

@author: nm2e16
"""


def positive_places(my_f, xs):
    """
    returns the positive numbers of the list. the returned numbers are positve
    after the function
    """
    res_list = []
    for num in xs:
        int_num = int(num)
        fun_num = my_f(int_num)
        if fun_num > 0:
            res_list.append(int_num)

    return res_list


def eval_f_0123(f):
    """
    returns the list of the functions after they have been through the function
    """
    res_list = []
    for i in range(4):
        value = f(i)
        print(value)
        res_list.append(value)
    return res_list
