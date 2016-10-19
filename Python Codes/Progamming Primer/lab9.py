# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 17:07:54 2016

@author: nm2e16
"""
import scipy.integrate
import math


def f(x):
    return x*x


def trapez(f, a, b, n):
    """
    returns the integral of a function using a trapeziodal message
    """
    h = (b-a)/n
    sum_int = 0
    for i in range(n-1):
        xi = a+((i+1)*h)
        sum_int = f(xi)+sum_int
    A = (h/2)*(f(a)+f(b)+(2*sum_int))
    return A


def finderror(n):
    A = trapez(f, 0, 2, n)
    return (8/3)-A


def using_quad():
    """returns the intergral and error using scipy"""
    return scipy.integrate.quad(f, 0, 2)


def std_dev(x):
    """returns the standard deviation of a list"""
    sum_xi = 0
    sigma_sum = 0
    for xi in x:
        sum_xi = sum_xi+xi

    mu = (1/len(x))*sum_xi
    print(mu)
    for xi in x:
        sigma_sum = (xi-mu)**2+sigma_sum
    sigma = math.sqrt(sigma_sum/len(x))
    return sigma


def reverse_dic(d):
    """returns a reversed dictionary"""
    new_values = [key for key in d.keys()]
    new_key = [val for val in d.values()]
    r = {}
    for i in range(len(new_key)):
        r[new_key[i]] = new_values[i]
    return r


def encode(code, msg):
    """
    returns a encrypted message
    """
    code_dict = code
    code_key = [let for let in code_dict.keys()]
    word = ''
    for let in msg:
        if let in code_key:
            let = code_dict[let]
        word = word+let
    return word


def code2():
    """A simple example i->s, s->g and g->i."""
    return {'i': 's', 's': 'g', 'g': 'i'}


def code3():
    """A more complicated code."""
    dic = {}
    # lower case letters
    dic['z'] = 'a'
    for i in range(ord('a'), ord('z')):
        dic[chr(i)] = chr(i + 1)
    # upper case letters
    dic['Z'] = 'A'
    for i in range(ord('A'), ord('Z')):
        dic[chr(i)] = chr(i + 1)
    # space, stop and some other special characters
    dic[' '] = '$'
    dic['.'] = '#'
    dic['#'] = '.'
    dic['$'] = ' '
    dic['?'] = '!'
    dic['!'] = '?'
    return dic


def decode(code, encoded_msg):
    """ returns a de-encrypted message"""
    word = ''
    new_code_enc = reverse_dic(code)
    code_key = [let for let in new_code_enc.keys()]

    for let in encoded_msg:
        if let in code_key:
            let = new_code_enc[let]
        word = word+let
    return word
