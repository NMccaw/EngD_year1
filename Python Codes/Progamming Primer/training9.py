"""Some support functions"""


def code0():
    """A trivial code - no change."""
    return {}


def code1():
    """A very simple example (symmetric)."""
    return {'e': 'x', 'x': 'e'}


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


def check_code_is_reversible(dic):
    """Given a dictionary used as a code mapping, this function checks
    whether the set of keys is the same set of values: if the elements
    in the keys are the same unique set as those in the values, then
    this mapping is bijective (the set of values is then actually a
    permutation of the set of input values) and can be inverted.

    If this is not the case, some debug information is printed, and a
    ValueError exception raised.
    """

    unique_keys = set()
    unique_vals = set()
    for key, val in dic.items():
        unique_keys.add(key)
        unique_vals.add(val)
    if unique_vals != unique_keys:
        print("Code is not reversible:")
        print("keys are   %s", sorted(unique_keys))
        print("values are %s", sorted(unique_vals))
        raise ValueError("Code is not reversible - stopping here")
    return True


def test_codes():
    """Check that codes defined above are reversible."""
    for c in (code0(), code1(), code2(), code3()):
        assert check_code_is_reversible(c)

def f(x):
    return x*x


def trapez(f,a,b,n):
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


def encode(code, msg):
    """
    returns a encrypted message
    """
    code_dict = code
    code_key = [let for let in code_dict.keys()]
    code_val = [let for let in code_dict.values()]
    word = ''
    for let in msg:
        if let in code_key:
            let = code_dict[let]
        word = word+let
    return word


def reverse_dic(d):
    """returns a reversed dictionary"""
    new_values = [key for key in d.keys()]
    new_key = [val for val in d.values()]
    r = {}
    for i in range(len(new_key)):
        r[new_key[i]] = new_values[i]
    return r
secretmessage = \
    "Zpv$ibwf$tvddfttgvmmz$efdpefe$uijt$tfdsfu$nfttbhf#$Dpohsbuvmbujpot?"


# if this file is executed on it's own, check codes given
if __name__ == "__main__":
    test_codes()
