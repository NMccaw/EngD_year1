# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 15:54:54 2016

@author: nm2e16
"""


def count_chars(s):
    """
    Returns the number of times characters occurs in a string in dictionary
    form
    """
    d = {}
    for char in s:
        if char not in d.keys():
            d[char] = s.count(char)
    return d



def derivative(f, x, eps=10**-6):
    """
    returns the central fininte difference of a function with variable accuracy
    """
    der = (f(x+(eps/2))-f(x-(eps/2)))/eps
    return der

def f(x):
     return x**2 -2


def g(x):
    return math.sin(x) +1.1


def newton(f, x, feps, maxit):
    """
    Returns the root of an equation via the newton method
    """
    count = 1
    while abs(f(x)) > feps:
        fprime = derivative(f, x, feps)
        x = x - f(x) / fprime
        count = count + 1
        if count == maxit:
            raise RuntimeError('Failed after %d iterations' % maxit)
    return x


def is_palindrome(s):
    """
    returns true or false if the string is a palindrome
    """
    if s == '':
        return True
    count = 1
    for i, char in enumerate(s):
        count = count + 1
        print(char)
    if char != s[-1-i]:
            return False
    return True
