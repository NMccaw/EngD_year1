import pytest
import class_test_2 as ct2

import numpy

def test_question_1():
    A_T = numpy.array([[3, 4], [2, 1]])
    assert(numpy.allclose(A_T, ct2.question_1()))

def test_question_2():
    v_squared = 253.5042735042735
    assert(numpy.allclose(v_squared, ct2.question_2()))

def test_question_3():
    w_even = 198.45327470248475
    assert(numpy.allclose(w_even, ct2.question_3()))

def test_question_4():
    x = numpy.array([0.18518519, -0.17283951,  0.74074074])
    assert(numpy.allclose(x, ct2.question_4()))

def test_question_5():
    minmax_evals = (0.423725877337, 206.027871927)
    assert(numpy.allclose(minmax_evals, ct2.question_5()))

def test_question_6():
    q_exact = numpy.array([[-0.87287156,  0.03539962, -0.48666426],
                     [-0.21821789, -0.92039002,  0.32444284],
                     [-0.43643578,  0.38939578,  0.81110711]])
    r_exact = numpy.array([[ -9.16515139,  -1.74574312,  -5.23722937],
                     [  0.        ,  13.45185418,  -4.99134589],
                     [  0.        ,   0.        ,   0.81110711]])
    q, r = ct2.question_6()
    assert(numpy.allclose(q_exact, q))
    assert(numpy.allclose(r_exact, r))

# Questions 7 and 8 will be checked visually

def test_question_9():
    root = 0.602449969894
    assert(numpy.allclose(root, ct2.question_9()))

def test_question_10():
    y5_exact = numpy.array([ 0.13606933,  0.47313955])
    y50_exact = numpy.array([ 0.11334685,  0.49624133])
    y5, y50 = ct2.question_10()
    assert(numpy.allclose(y5_exact, y5))
    assert(numpy.allclose(y50_exact, y50))

if __name__ == "__main__":
    pytest.main()
    ct2.question_7()
    ct2.question_8()
