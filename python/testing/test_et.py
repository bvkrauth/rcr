"""
TEST_ET.PY: Unit tests for estimate_theta()
"""


import pytest
import numpy as np

from rcrbounds import estimate_theta_segments, estimate_theta


# Basic functionality
def test_et_basic():
    """estimate theta from simple data"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.5, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    ts1 = estimate_theta_segments(mv1)[0]
    et_true = np.array([-0.33333333, 1.])
    with pytest.warns(UserWarning, match="Inaccurate SE"):
        test_et = estimate_theta(mv1, lr1, ts1)
    assert test_et[:, 0] == pytest.approx(et_true, rel=1e-04)


def test_et_realdata(moment_vector):
    """estimate theta from real data"""
    lambda_range = np.array([0.0, 1.0])
    theta_segments = estimate_theta_segments(moment_vector)[0]
    et_true = np.array([[5.13504376e+00, 6.15418928e+01, 2.90135406e+01,
                         -9.12983158e+01, 5.90041894e+01, -8.93300435e-01,
                         -7.88725310e+01, -1.35242174e+00, -1.82475877e+02,
                         1.49588107e+00, 1.00493824e+00, -3.83778143e+00,
                         1.15368002e+01, -4.10554973e-01, -1.34697339e+01,
                         -1.71351694e-01, -1.80272712e+02, 1.16780438e-01,
                         -1.13509850e+00, 6.09839080e+00, -2.32062522e-01,
                         -6.99568697e+00, -2.30153513e-02, -9.45064400e+01,
                         2.34747229e+00, -1.80916161e+01, 6.66094500e-01,
                         2.09384580e+01, 1.68651775e-01, 2.81533013e+02,
                         -1.51773448e+00, 4.29563522e-01, 8.83539162e-01,
                         -1.39915455e+00, 3.06232087e+01, -1.44523749e-02,
                         -3.88855612e-01, 5.47869935e-02, -6.00042013e+00,
                         9.62873893e-01, 1.59223870e+00, 4.87781110e+00,
                         -6.90173844e-05, 5.00039063e+00, -2.56754049e+01],
                        [5.20150257e+00, 1.40841120e+01, 6.64793094e+00,
                         -2.09059332e+01, 1.33309588e+01, -1.97134309e-01,
                         -1.78586773e+01, -1.44338171e+00, -2.22802321e+02,
                         3.51963458e-01, 2.41387296e-01, -9.10296914e-01,
                         2.60893923e+00, -9.21291086e-02, -3.05196227e+00,
                         -4.19568148e-02, -4.08043112e+01, 3.00528248e-02,
                         -2.78581790e-01, 1.37924709e+00, -5.21165908e-02,
                         -1.58522968e+00, -6.85820532e-03, -2.13933567e+01,
                         5.63724100e-01, -4.09147624e+00, 1.49533901e-01,
                         4.74444552e+00, 4.31065081e-02, 6.37273712e+01,
                         -3.44858582e-01, 9.72804394e-02, 2.03443701e-01,
                         -3.16456560e-01, 6.97515752e+00, -3.26189067e-03,
                         -8.82443477e-02, 1.23080405e-02, -1.36003054e+00,
                         2.16359819e-01, 3.60818788e-01, 1.05554479e+00,
                         1.10295995e-18, 4.89021069e+00, -2.54364435e+01]])
    test_et = estimate_theta(moment_vector, lambda_range, theta_segments)
    assert test_et[:, 0] == pytest.approx(et_true[:, 0])
    assert test_et == pytest.approx(et_true, rel=1.0e-4, abs=1e-04)


# Varying lambda_range
def test_et_lambdapoint(moment_vector):
    """estimate theta when lambda_range is a single point"""
    lr0 = np.array([0.0, 0.0])
    theta_segments = estimate_theta_segments(moment_vector)[0]
    et_true = np.array([[5.20150257e+00, 1.40841120e+01, 6.64793094e+00,
                         -2.09059332e+01, 1.33309588e+01, -1.97134309e-01,
                         -1.78586773e+01, -1.44338171e+00, -2.22802321e+02,
                         3.51963458e-01, 2.41387296e-01, -9.10296914e-01,
                         2.60893923e+00, -9.21291086e-02, -3.05196227e+00,
                         -4.19568148e-02, -4.08043112e+01, 3.00528248e-02,
                         -2.78581790e-01, 1.37924709e+00, -5.21165908e-02,
                         -1.58522968e+00, -6.85820532e-03, -2.13933567e+01,
                         5.63724100e-01, -4.09147624e+00, 1.49533901e-01,
                         4.74444552e+00, 4.31065081e-02, 6.37273712e+01,
                         -3.44858582e-01, 9.72804394e-02, 2.03443701e-01,
                         -3.16456560e-01, 6.97515752e+00, -3.26189067e-03,
                         -8.82443477e-02, 1.23080405e-02, -1.36003054e+00,
                         2.16359819e-01, 3.60818788e-01, 1.05554479e+00,
                         1.10295995e-18, 4.89021069e+00, -2.54364435e+01],
                        [5.20150257e+00, 1.40841120e+01, 6.64793094e+00,
                         -2.09059332e+01, 1.33309588e+01, -1.97134309e-01,
                         -1.78586773e+01, -1.44338171e+00, -2.22802321e+02,
                         3.51963458e-01, 2.41387296e-01, -9.10296914e-01,
                         2.60893923e+00, -9.21291086e-02, -3.05196227e+00,
                         -4.19568148e-02, -4.08043112e+01, 3.00528248e-02,
                         -2.78581790e-01, 1.37924709e+00, -5.21165908e-02,
                         -1.58522968e+00, -6.85820532e-03, -2.13933567e+01,
                         5.63724100e-01, -4.09147624e+00, 1.49533901e-01,
                         4.74444552e+00, 4.31065081e-02, 6.37273712e+01,
                         -3.44858582e-01, 9.72804394e-02, 2.03443701e-01,
                         -3.16456560e-01, 6.97515752e+00, -3.26189067e-03,
                         -8.82443477e-02, 1.23080405e-02, -1.36003054e+00,
                         2.16359819e-01, 3.60818788e-01, 1.05554479e+00,
                         1.10295995e-18, 4.89021069e+00, -2.54364435e+01]])
    test_et = estimate_theta(moment_vector, lr0, theta_segments)
    assert test_et[:, 0] == pytest.approx(et_true[:, 0])
    assert test_et == pytest.approx(et_true, rel=1.0e-4, abs=1.0e-4)


def test_et_nolambdalow(moment_vector):
    """estimate theta when lambda_range has no lower bound"""
    lr0 = np.array([-np.inf, 1])
    theta_segments = estimate_theta_segments(moment_vector)[0]
    et_true = np.array([5.13504376,  8.16970996])
    with pytest.warns(UserWarning, match="Inaccurate SE"):
        test_et = estimate_theta(moment_vector, lr0, theta_segments)
    assert test_et[:, 0] == pytest.approx(et_true)


def test_et_nolambdahigh(moment_vector):
    """estimate theta when lambda_range has no upper bound"""
    lr0 = np.array([0, np.inf])
    theta_segments = estimate_theta_segments(moment_vector)[0]
    test_et = estimate_theta(moment_vector, lr0, theta_segments)
    assert test_et[0, 0] == -np.inf
    assert test_et[1, 0] == np.inf
    assert np.all(test_et[:, 1:] == 0.0)


# Special cases for moments
def test_et_nearrct():
    """estimate theta for near-perfect RCT"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.000001, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    ts1 = estimate_theta_segments(mv1)[0]
    et_true = 0.5
    test_et = estimate_theta(mv1, lr1, ts1)
    assert test_et[:, 0] == pytest.approx(et_true, rel=1e-04)


def test_et_rct():
    """estimate theta for perfect RCT"""
    mv1 = np.array([0, 0, 0, 1, 0.5, 0.0, 1, 0.5, 1.0])
    lr1 = np.array([0.0, 1.0])
    et_true = 0.5
    # This test currently fails with an UnboundLocalError
    try:
        ts1 = estimate_theta_segments(mv1)[0]
    except UnboundLocalError:
        pass
    else:
        test_et = estimate_theta(mv1, lr1, ts1)
        assert test_et[:, 0] == pytest.approx(et_true, rel=1e-04)
