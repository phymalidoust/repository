from ase.clease import LinearRegression
import numpy as np
from scipy.special import polygamma
from scipy.optimize import brentq
import time
import json
from itertools import product

class BayesianCompressiveSensing(LinearRegression):
    """
    Fit a sparse CE model to data. Based on the method 
    described in

    Babacan, S. Derin, Rafael Molina, and Aggelos K. Katsaggelos. 
    "Bayesian compressive sensing using Laplace priors." 
    IEEE Transactions on Image Processing 19.1 (2010): 53-63.

    Different values has different priors.
    1. For the ECIs a normal distribution is assumed
        (the i-th eci is: eci_i -- N(J | 0, var_i)=
    2. The inverce variance of each ECI is gamma distributed
        (i.e. 1/var_i -- gamma(x | 1, lambda/2))
    3. The lambda parameter above is also gamma distributed
        (i.e. lamb -- gamma(x | shape_lamb/2, shape_lamb/2))
    4. The noise parameter is uniformly distributed on the 
        positive axis (i.e. noise -- uniform(x | 0, inf)

    Arguments
    ===========
    shape_var: float
        Shape parameter for the gamma distribution for the 
        inverse variance (1/var -- gamma(x | shape_var/2, rate_var/2))
    rate_var: float
        Rate parameter for the gamma distribution for the 
        inverse variance (1/var -- gamma(x | shape_var/2, rate_var/2))
    shape_lamb: float
        Shape parameter for gamma distribution for the 
        lambda parameter (lambda -- gamma(x | 1, shape_lamb))
    variance_opt_start: int
        Optimization of lambda, shape_lamb and 
        inverse variance starts after this amount of iterations
    fname: str
        Backup file for parameters
    maxiter: int
        Maximum number of iterations
    output_rate_sec: int
        Interval in seconds between status messages
    select_strategy: str
        Strategy for selecting new correlation function for each iteration. 
        If 'max_increase' it will select the basis function that leads
        to the largest increase in likelihood value. 
        If 'random' correlation functions are selected at random
    noise: float
        Initial estimate of the noise in the data
    """
    def __init__(self, shape_var=0.5, rate_var=0.5, shape_lamb=0.5, 
                 variance_opt_start=100, fname="bayes_compr_sens.json",
                 maxiter=100000, output_rate_sec=2,
                 select_strategy="max_increase", noise=0.1):
        LinearRegression.__init__(self)

        # Paramters
        self.shape_var = shape_var
        self.rate_var = rate_var
        self.shape_lamb = shape_lamb
        self.variance_opt_start = variance_opt_start
        self.maxiter = maxiter
        self.output_rate_sec = output_rate_sec
        self.select_strategy = select_strategy
        self.fname = fname
        self.noise = noise

        # Arrays used during fitting
        self.X = None
        self.y = None

        self.gammas = None
        self.eci = None
        self.inv_variance = None
        self.lamb = None
        self.inverse_sigma = None
        self.eci = None

        # Quantities used for fast updates
        self.S = None
        self.Q = None
        self.ss = None
        self.qq = None

    def _initialize(self):
        """
        Initialize all parameters after X and y is given
        """
        num_features = self.X.shape[1]
        if self.gammas is None:
            self.gammas = np.zeros(num_features)
        self.eci = np.zeros_like(self.gammas)

        if self.inv_variance is None:
            self.inv_variance = 1.0/self.noise**2
        
        if self.lamb is None:
            self.lamb = 0.0
            #self.lamb = self.optimal_lamb()

        self.inverse_sigma = np.zeros((num_features, num_features))
        self.eci = np.zeros(num_features)

        if self.X is not None:
            # Quantities used for fast updates
            self.S = np.diag(self.inv_variance*self.X.T.dot(self.X))
            self.Q = self.inv_variance*self.X.T.dot(self.y)
            self.ss = self.S/(1.0 - self.gammas*self.S)
            self.qq = self.Q/(1.0 - self.gammas*self.S)

    def precision_matrix(self, X):
        """
        Return the precision matrix needed by the Evaluate class.
        Only contributions from the correlation functions with 
        gamma > 0 are included.
        """
        if not np.allclose(X, self.X):
            raise RuntimeError("Inconsistent design matrix given!")
        sel_indx = self.selected
        X_sel = self.X[:, sel_indx]
        prec = np.linalg.inv(X_sel.T.dot(X_sel))

        N = self.X.shape[1]
        full_prec = np.zeros((N, N))

        indx = list(range(len(sel_indx)))
        for i in product(indx, indx):
            full_prec[sel_indx[i[0]], sel_indx[i[1]]] = prec[i[0], i[1]]
        return full_prec

    def mu(self):
        """
        Calculate the expectation value for the ECIs
        """
        sel = self.selected
        return self.inv_variance*self.inverse_sigma.dot(self.X[:, sel].T.dot(self.y))

    def optimal_gamma(self, indx):
        """
        Return the gamma value that maximize the likelihood

        Arguments
        ==========
        indx: int
            Index of the selected correlation function
        """
        s = self.ss[indx]
        qsq = self.qq[indx]**2

        if self.lamb < 1E-6:
            return (qsq - s)/s**2
        term1 = s + 2*self.lamb

        delta = s**2 + 4*self.lamb*qsq

        gamma = (np.sqrt(delta) - term1)/(2*self.lamb*s)
        assert np.sign(gamma) == np.sign(qsq - s - self.lamb)
        return gamma

    def optimal_lamb(self):
        """Calculate the optimal value for the lambda parameter. """
        N = self.X.shape[1]
        return 2*(N - 1 + 0.5*self.shape_lamb)/(np.sum(self.gammas) + self.shape_lamb)

    def optimal_inv_variance(self):
        """Calculate the optimal value for the inverse variance"""
        N = self.X.shape[1]
        a = 1.0
        b = 0.0
        mse = np.sum((self.y - self.X.dot(self.eci))**2)
        return (0.5*N + a)/(0.5*mse + b)

    def optimal_shape_lamb(self):
        """Calculate the optimal value for the shape paremeter for lambda."""
        res = brentq(shape_parameter_equation, 1E-30, 1E100, args=(self.lamb,), maxiter=10000)
        return res

    def update_quantities(self):
        """Update helper parameters needed for the next iteration."""
        sel = self.selected
        X_sel = self.X[:, sel]
        prec = X_sel.dot(self.inverse_sigma).dot(X_sel.T)


        self.S = np.diag(self.inv_variance*self.X.T.dot(self.X) - self.inv_variance**2 * self.X.T.dot(prec.dot(self.X)))
        self.Q = self.inv_variance*self.X.T.dot(self.y) - self.inv_variance**2 * self.X.T.dot(prec.dot(self.y))

        self.ss = self.S/(1.0 - self.gammas*self.S)
        self.qq = self.Q/(1.0 - self.gammas*self.S)

    @property
    def selected(self):
        return np.argwhere(self.gammas > 0.0)[:, 0]

    def update_sigma_mu(self):
        """Update sigma and mu."""
        X_sel = self.X[:, self.selected]
        self.inverse_sigma = np.linalg.inv(self.inv_variance*X_sel.T.dot(X_sel) + np.diag(1.0/self.gammas[self.selected]))
        self.eci[self.selected] = self.mu()

    def get_basis_function_index(self, select_strategy):
        """Select a new correlation function."""
        if select_strategy == "random":
            return np.random.randint(low=0, high=len(self.gammas))
        elif select_strategy == "max_increase":
            return self._get_bf_with_max_increase()

    def log_likelihood_for_each_gamma(self, gammas):
        """Log likelihood value for all gammas.

        Arguments
        =========
        gammas: np.ndarray
            Value for all the gammas
        """
        denum = 1+gammas*self.ss
        if np.any(denum < 0.0):
            raise RuntimeError("Negative denuminator inside log function. Try to "
                               "increase the noise parameter.")
        return np.log(1/denum) + self.qq**2*gammas/denum - self.lamb*gammas

    def _get_bf_with_max_increase(self):
        """Return the index of the correlation function that leads
            to the largest increase in likelihiood value."""

        new_gammas = np.array([self.optimal_gamma(i) for i in range(len(self.gammas))])
        if np.all(new_gammas < 0.0):
            raise RuntimeError("All gammas are smaller than 0! Cannot include any!")

        new_gammas[new_gammas < 0.0] = 0.0
        current_likeli = self.log_likelihood_for_each_gamma(self.gammas)
        new_likeli = self.log_likelihood_for_each_gamma(new_gammas)

        diff = new_likeli - current_likeli
        return np.argmax(diff)

    def rmse(self):
        """Return root mean square error."""
        indx = self.selected
        pred = self.X[:, indx].dot(self.eci[indx])
        return np.sqrt(np.mean((pred - self.y)**2))

    def log(self, msg):
        print(msg)

    @property
    def num_ecis(self):
        return np.count_nonzero(self.gammas)

    def to_dict(self):
        """Convert all parameters to a dictionary."""
        data = {}
        data["inv_variance"] = self.inv_variance
        data["gammas"] = self.gammas.tolist()
        data["shape_var"] = self.shape_var
        data["rate_var"] = self.rate_var
        data["shape_lamb"] = self.shape_lamb
        data["lamb"] = self.lamb
        data["maxiter"] = self.maxiter
        data["output_rate_sec"] = self.output_rate_sec
        data["select_strategy"] = self.select_strategy
        data["noise"] = self.noise
        return data

    def save(self):
        """Save the results from file."""
        with open(self.fname, 'w') as outfile:
            json.dump(self.to_dict(), outfile)
        print("Backup data written to {}".format(self.fname))

    @staticmethod
    def load(fname):
        bayes = BayesianCompressiveSensing()
        bayes.fname = fname
        with open(fname, 'r') as infile:
            data = json.load(infile)
        
        bayes.inv_variance = data["inv_variance"]
        bayes.gammas = np.array(data["gammas"])
        bayes.shape_var = data["shape_var"]
        bayes.rate_var = data["rate_var"]
        bayes.shape_lamb = data["shape_lamb"]
        bayes.lamb = data["lamb"]
        bayes.maxiter = data["maxiter"]
        bayes.output_rate_sec = data["output_rate_sec"]
        bayes.select_strategy = data["select_strategy"]
        bayes.noise = data["noise"]
        return bayes
    
    def __eq__(self, other):
        """Compare to BayesianCompressiveSensing objects."""
        equal = True

        # Required fields to be equal if two objects
        # should be considered equal
        items = ["fname", "gammas", "inv_variance", "lamb",
                 "shape_var", "rate_var", "shape_lamb", "lamb",
                 "maxiter", "select_strategy", "output_rate_sec",
                 "noise"]
        for k in items:
            v = self.__dict__[k]
            if isinstance(v, np.ndarray):
                equal = equal and np.allclose(v, other.__dict__[k])
            elif isinstance(v, float):
                equal = equal and abs(v - other.__dict__[k]) < 1E-6
            else:
                equal = equal and (v == other.__dict__[k])
        return equal

    def estimate_loocv(self):
        """Return an estimate of the LOOCV."""
        X_sel = self.X[:, self.selected]
        e_pred = self.X.dot(self.eci)
        delta_e = e_pred - self.y
        prec = np.linalg.inv(X_sel.T.dot(X_sel))
        cv_sq = np.mean((delta_e / (1 - np.diag(X_sel.dot(prec).dot(X_sel.T))))**2)
        return np.sqrt(cv_sq)

    def fit(self, X, y):
        """Fit ECIs to the data

        Arguments
        =========
        X: np.ndarray
            Design matrix (NxM: N number of datapoints, M number of correlation functions)
        y: np.ndarray
            Array of length N with the energies
        """
        allowed_strategies = ["random", "max_increase"]

        if self.select_strategy not in allowed_strategies:
            raise ValueError("select_strategy has to be one of {}"
                             "".format(allowed_strategies))

        self.X = X
        self.y = y
        self._initialize()

        is_first = True
        iteration = 0
        now = time.time()
        d_gamma = 1E100
        while iteration < self.maxiter:
            if time.time() - now > self.output_rate_sec:
                msg = "Iter: {} ".format(iteration)
                msg += "RMSE: {:.3E} ".format(1000.0*self.rmse())
                msg += "LOOCV (approx.): {:.3E}".format(1000.0*self.estimate_loocv())
                msg += "Num ECI: {} ".format(self.num_ecis)
                msg += "Lamb: {:.3E}. Shape lamb: {:.3E} ".format(self.lamb, self.shape_lamb)
                msg += "Noise: {:.3E}".format(np.sqrt(1.0/self.inv_variance))
                self.log(msg)
                now = time.time()

            iteration += 1
            already_excluded = False

            if is_first:
                indx = np.argmax(self.qq**2 - self.ss)
                is_first = False
            else:
                indx = self.get_basis_function_index(self.select_strategy)

            gamma = self.optimal_gamma(indx)
            d_gamma = gamma - self.gammas[indx]
            if gamma > 0.0:
                self.gammas[indx] = gamma
            else:
                gamma = self.gammas[indx]

                if abs(gamma) < 1E-6:
                    already_excluded = True

                self.gammas[indx] = 0.0
                self.eci[indx] = 0.0

            if already_excluded:
                continue
            
            self.update_sigma_mu()
            self.update_quantities()

            if iteration > self.variance_opt_start:
                self.lamb = self.optimal_lamb()
                self.shape_lamb = self.optimal_shape_lamb()
                self.inv_variance = self.optimal_inv_variance()

            if abs(d_gamma) < 1E-8:
                break

        # Save backup for future restart
        if self.fname:
            self.save()
        return self.eci

    def show_shape_parameter(self):
        """Show a plot of the transient equation for the optimal
            shape parameter for lambda."""
        from matplotlib import pyplot as plt
        x = np.logspace(-10, 10)

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(x, shape_parameter_equation(x, self.lamb))
        ax.axhline(0, ls="--")
        ax.set_xscale("log")
        plt.show()


def shape_parameter_equation(x, lamb):
    return np.log(x/2.0) + 1 - polygamma(0, x/2) + np.log(lamb) - lamb