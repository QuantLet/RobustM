import numpy as np
from math import pi, sqrt, exp, log
from scipy.sparse.linalg import eigsh

from cvxopt import solvers, matrix, spmatrix
solvers.options['show_progress'] = False

def block_iterator(n, block_num):
    m = n // block_num
    r = n % block_num
    
    lst = []
    for j in range(block_num):
        if j < r:
            a = j * (m + 1)
            b = (j + 1) * (m + 1)
        else:
            a = r * (m + 1) + (j - r) * m
            b = r * (m + 1) + (j - r + 1) * m
        lst.append((a, b))
    
    return lst

def project(w):
    d = np.size(w)
    return w + np.ones(d) * (1 - np.sum(w)) / d

def project_bregman(x, delta):
    n = np.size(x)
    y = np.sort(x)
    s = 0
    C = 0
    for j, (u, u_next) in enumerate(zip(y[:-1], y[1:])):
        if u <= 0:
            continue
        
        s += u
        C_ = (1 - delta * (n - j)) / s
        if C_ < 0:
            continue
        if C_ * u_next >= delta:
            C = C_
            break
    if C == 0:
        C = 1 / (s + u_next)
    return np.maximum(np.minimum(x * C, delta), 0)

def top_eigenvector(M):
    n, _ = np.shape(M)
    vals, vecs = eigsh(M, k=1)
    return vals[0], vecs[:, 0]
    

def hopkins_mean(x, block_num, eta=0.1, steps=50):
    n, d = np.shape(x)
    z = np.array([np.mean(x[a:b, :], axis=0) for (a, b) in block_iterator(n, block_num)])
    
    return hopkins_mean_from_buckets(z, eta=eta, steps=steps)

def hopkins_mean_from_buckets(z, eta=0.1, steps=50):
    block_num, d = np.shape(z)
    
    vbest = None
    normbest = -1
    
    w = np.ones(block_num) / block_num
    for _ in range(steps):
        v = np.dot(w, z)
        diff = z - v
        M = np.matmul(np.matmul(diff.T, np.diag(w)), diff)
        val, vec = top_eigenvector(M)
        if val < normbest or normbest == -1:
            normbest = val
            vbest = v
        
        tau = np.square(np.dot(diff, vec))
        eta_careful = min(eta, 1 / (2 * np.max(tau)))
        w = w * (np.ones(block_num) - eta_careful * tau)
        w = project_bregman(w, 4 / (3 * block_num))
    return vbest

def robust_qp(X, block_num, lmbd, gamma, w0=None, steps=100, project_steps=None):
    n, d = np.shape(X)
    
    covs = np.array([
        np.cov(X[a:b, :].T) for (a, b) in block_iterator(n, block_num)
    ])
    mu = hopkins_mean(X, block_num)
    
    if w0 is None:
        #k = min(2, d // 2)
        #w0 = np.append(np.zeros(d - k), np.ones(k) / k)
        w0 = np.ones(d) / d
    
    w = np.array(w0)
    for _ in range(steps):
        w = w - gamma * (hopkins_mean_from_buckets(np.dot(covs, w)) - lmbd * mu)
        w = project(w)
    
    return w

def recommended_gamma(emp_cov):
    return 1 / (4 * np.linalg.norm(emp_cov, ord=2))

def plain_qp(mu, sigma, lmbd):
    d = np.size(mu)
    ones = np.ones(d)
    mu_inv = np.linalg.solve(sigma, mu)
    ones_inv = np.linalg.solve(sigma, ones)
    
    delta = (1 - lmbd * np.dot(ones, mu_inv)) / np.dot(ones, ones_inv)
    return lmbd * mu_inv + delta * ones_inv

def plain_ge_qp(mu, sigma, lmbd, ge):
    assert(ge >= 1.0)
    d = np.size(mu)
    
    eye = spmatrix(1.0, range(d), range(d))
    Omat = spmatrix(0.0, range(d), range(d))
    ones_d = matrix(1.0, (1, d))
    zeros_d = matrix(0.0, (1, d))
    
    P = matrix([[matrix(sigma), Omat], [Omat, Omat]])
    q = matrix([matrix(-lmbd * mu), zeros_d.T])
    
    G = matrix([[eye, -eye, zeros_d], [-eye, -eye, ones_d]])
    h = matrix([zeros_d.T, zeros_d.T, ge])
    
    A = matrix([[ones_d], [zeros_d]])
    b = matrix(1.0, (1,1))
    
    sol = solvers.qp(P, q, G, h, A, b)
    return np.reshape(np.array(sol['x'][:d]), (d,))
