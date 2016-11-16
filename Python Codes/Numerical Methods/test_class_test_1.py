import pytest
import class_test_1 as ct1

import numpy

def test_question_1():
    A_cubed = numpy.array([[85.0, 129.0], [172.0, 257.0]])
    assert(numpy.allclose(A_cubed, ct1.question_1()))

def test_question_2():
    v_squared = 93.504273504273499
    assert(numpy.allclose(v_squared, ct1.question_2()))

def test_question_3():
    w_even = 11.354570156340603
    assert(numpy.allclose(w_even, ct1.question_3()))

def test_question_4():
    x = numpy.array([-1.26, 0.64, 3.2])
    assert(numpy.allclose(x, ct1.question_4()))

def test_question_5():
    minmax_evals = (1.3490577004552573, 16.353671026158455)
    assert(numpy.allclose(minmax_evals, ct1.question_5()))

def test_question_6():
    A_L = numpy.array([[0, 0, 0], [2, 0, 0], [4, 6, 0]])
    assert(numpy.allclose(A_L, ct1.question_6()))

# Questions 7 and 8 will be checked visually

def test_question_9():
    root = 0.47587757384907897
    assert(numpy.allclose(root, ct1.question_9()))

def test_question_10():
    y5_exact = numpy.array([ 0.67882719,  0.70681336])
    y50_exact = numpy.array([ 0.76816916,  0.69481969])
    y5, y50 = ct1.question_10()
    assert(numpy.allclose(y5_exact, y5))
    assert(numpy.allclose(y50_exact, y50))

if __name__ == "__main__":
    pytest.main()
    ct1.question_7()
    ct1.question_8()
