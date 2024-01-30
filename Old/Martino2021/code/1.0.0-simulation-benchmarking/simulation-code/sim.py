import unittest
import numpy as np
from scipy.signal import medfilt
from scipy.special import softmax
from skbio.stats.composition import closure
from numpy.random import poisson, lognormal, normal, randint


def softmax_inv(x):
    return np.log(x)


class TestSoftInv(unittest.TestCase):

    def setUp(self):
        self.exp = closure(np.array([[100, 5, 20],
                                     [200, 234, 14]]))
        pass

    def testinv(self):

        res = softmax(softmax_inv(self.exp))
        res *= 1 / res.sum(axis=1)[0]
        try:
            np.testing.assert_array_almost_equal(res, self.exp, 6)
        except AssertionError:
            print('arrays not equal')


def simulate(X, depth=100, read_std=10.0,
             pnormal=0.1, pheter=0.1,
             pzero=0.1, mwindow=False):
    """
    Poisson Log-Normal simulation
    one real count data. The proceedure
    follows the same steps as
    Äijö et al. [1] but with a Poisson
    Log-Normal count model. This allows
    the total read depth to vary across
    samples which can be modeled by the
    dispersion perameter of the Poisson.

    Parameters
    ----------
    table : array like
        A matrix of counts.
    depth : int
        Mean read depth of the
        simulation.
    read_std : float
        Controls dispersion.
    pnormal : float [0,1]
        Percent
        normally-distributed
        noise to the data.
    pheter : float [0,1]
        Percent random
        noise to the data.
    pzero : float [0,1]
        Percent zeros added.

    Returns
    -------
    sim : array-like
        simulated counts.

    References
    ----------
    .. [1] Äijö, Tarmo, Christian L. Müller,
           and Richard Bonneau. 2018.
           “Temporal Probabilistic Modeling of Bacterial
           Compositions Derived from 16S rRNA Sequencing.”
           Bioinformatics  34 (3): 372–80.
    """

    # subset and order single digester
    sim = X.copy()
    # (1) estimate proportions
    sim = closure(sim)
    # (2) filter low-mean reads
    abund_filt_mask = sim.mean(axis=0) > 1e-5
    mask_index = [i for i, v in enumerate(abund_filt_mask) if v]
    sim = sim.T[abund_filt_mask].T
    # (3) median filter by feature
    if mwindow:
        for i in range(sim.shape[1]):
            sim[:, i] = medfilt(sim[:, i], mwindow)
    sim = closure(sim)  # check closure
    # (4) project onto real
    # with iverse softmax
    sim = softmax_inv(sim)
    # (5) noise gauss (@10% of samps.)
    err = pnormal * np.ones_like(sim)
    sim = normal(sim, err)
    # add some random noise (@10% of samps.)
    err = pheter * np.ones_like(sim)
    i = randint(0, err.shape[0], 5000)
    j = randint(0, err.shape[1], 5000)
    err[i, j] = 0.5
    sim = normal(sim, err)
    # add random zeros (@10% of samps.)
    err = pzero * np.ones_like(sim)
    i = randint(0, err.shape[0], 100)
    j = randint(0, err.shape[1], 100)
    err[i, j] = 0.0
    sim = normal(sim, err)
    # (6) back to prop
    sim = closure(softmax(sim))
    # (7) Poisson - Log-Normal Counts
    mean = depth
    stdev = read_std
    phi = (stdev ** 2 + mean ** 2) ** 0.5
    mu = mean**2 / phi
    mu = np.log(mu * sim.T)
    sigma = (np.log(phi ** 2 / mean ** 2)) ** 0.5
    sim = np.vstack([poisson(lognormal(mu[:, i],
                                       sigma))
                     for i in range(sim.shape[0])])

    return sim, mask_index, sigma
