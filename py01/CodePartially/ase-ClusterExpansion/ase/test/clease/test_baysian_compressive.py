import os
from ase.clease import BayesianCompressiveSensing
import numpy as np
from scipy.special import polygamma

# Construct 30 measurements with 400 basis functions
np.random.seed(0)
fname = "test_bayes_compr_sens.json"

bayes = BayesianCompressiveSensing(fname=fname, output_rate_sec=2, 
                                   maxiter=100)
   
def test_optimize_shape_parameter(bayes):
    bayes.lamb = 1.0
    opt = bayes.optimal_shape_lamb()
    assert abs(np.log(opt/2.0) - polygamma(0, opt/2.0)) < 1E-6

def test_fit(bayes):
    X = np.random.rand(30, 400)
    y = 60.0*X[:, 20] - 80.0*X[:, 2]
    bayes.fit(X, y)

    expected_eci = np.zeros(X.shape[1])
    expected_eci[20] = 60.0
    expected_eci[2] = -80.0
    assert np.allclose(bayes.eci, expected_eci, rtol=1E-4)

def test_fit_more_coeff():
    bayes = BayesianCompressiveSensing(fname=fname, noise=0.1)
    X = np.random.rand(30, 400)
    coeff = [6.0, -2.0, 5.0, 50.0, -30.0]
    indx = [0, 23, 19, 18, 11]
    y = 0.0
    expected_eci = np.zeros(X.shape[1])
    for c, i in zip(coeff, indx):
        y += X[:, i]*c
        expected_eci[i] = c
    bayes.fit(X, y)
    assert np.allclose(bayes.eci, expected_eci, atol=1E-2)
    

def test_save_load(bayes):
    bayes.save()

    bayes2 = BayesianCompressiveSensing.load(fname)
    assert bayes == bayes2

test_optimize_shape_parameter(bayes)
test_fit(bayes)
test_save_load(bayes)
test_fit_more_coeff()
os.remove(fname)
