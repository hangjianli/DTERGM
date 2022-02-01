import numpy as np
import time
import logging

logging.basicConfig(level=logging.INFO)

def gflasso(y, lams, verbose=0):
    """"""
    n, p = y.shape
    
    mean_signal = np.mean(y, axis=0)
    w = default_weights(n)
    lams = sorted(lams, key=lambda x: -x)
    nlambda = len(lams)
    C = left_multiply_by_xt(y, w)
    
    beta = []
    active_set = []
    res = {}
    for lam in lams:
        # Optimize for this lambda with a warm restart from the previous lambda
        output = optimize(beta, active_set, lam, C, y, w, verbose=verbose)
        res[lam] = output.copy()
        
    return res

def onetwo_norm(A):
    return np.sqrt(np.sum(A**2, axis=0)).sum()

def default_weights(n):
    """Generate default weights according to Eq.(5) in Bleakley & Vert, 2011 """
    a = 1 + np.arange(n-1)
    w = np.sqrt(n / (a * (n-a)))
    return w


def left_multiply_by_xt(Y, w):
    """Fast computation for X'*Y for the group fused Lasso"""
    n, p = Y.shape
    u = np.cumsum(Y, axis=0)
    part1 = np.arange(1,n).reshape(-1,1) @ u[-1].reshape(1,-1) / n - u[:-1]
    part2 = np.array([w,]*p).transpose()
    C = np.multiply(part1, part2)
    return C


def reconstruct_signal(n, act_set, beta, w, meansignal):
    """Reconstruct a piecewise constant profile from its increments"""
    a, p = beta.shape 
    # a: number of change points
    # p: dimension of signal
    
    signal = np.zeros((n, p))
    weights = np.array([w[act_set.reshape(-1,1) - 1]] * p).transpose()
    signal[act_set.flatten()] = np.multiply(beta, weights).reshape(-1, p)
    signal = np.cumsum(signal, axis=0)
    
    # set the mean signla to meansignal
    m = meansignal - np.mean(signal, axis=0)
    signal += np.array([m]*n)
    
    return signal


def xtx(n, A, B, w):
    """
    Compute X[:, A].T * X[:, B], where X is the design matrix of the weighted GFL
    
    Parameters:
    ----------
    n:
    A: a column vector of indices in [0, n-2]
    B: a column vector of indices in [0, n-2]
    w: (n-1) * 1 column vector of weights
    """
    
    assert type(A) is np.ndarray, f"arg A: {A} is not np.ndarray!"
    assert type(B) is np.ndarray, f"arg B: {B} is not np.ndarray!"
    
    u = np.array([A+1]*len(B)).transpose().reshape(len(A), len(B))
    v = np.array([B+1]*len(A)).reshape(len(A), len(B))
    tmp = np.multiply(n-np.maximum(u,v), w[A] @ w[B].T)/n
    g = np.multiply(np.minimum(u, v), tmp)
    return g
    
    
def block_coordinate_descent(beta, active_set, xty, n, w, lam, verbose=True):
    """Block coordinate descent over the active set of beta
    
    Parameters
    ----------
    beta: (np.array): a by p matrix as initializer
    active_set: (array): a by 1 matrix with indices of nonzero beta_t, in the range of [1, T-1]
    
    
    Returns:
            
    """
    maxit = 1e5
    
    if len(beta) == 0:
        return beta
    
    a, p = beta.shape
    tol = 1e-8 * p
    norm_beta = np.linalg.norm(beta, axis=1)
    gain = 2*tol*np.ones(a);
    
    itc = 0
    while(any(gain > tol) & (itc < maxit)):
        i = np.mod(itc, a)
        asi = active_set[i]

        # Compute the vector S
        XitX = xtx(n, asi, active_set, w)
        gammai = XitX.flatten()[i]
        # the indices of the active set excluding the one being optimized. This could be empty
        indwithouti = np.array([k for k in range(a) if k != i])
        if len(indwithouti) == 0:
            S = xty[i, :]
        else:
            S = xty[i, :] - XitX.flatten()[indwithouti] @ beta[indwithouti, :]
        nS = np.linalg.norm(S)
        if nS < lam:
            newbeta = np.zeros(p) # If |S| < lambda, then softshriholding returns zero
        else:
            newbeta = (S * (1 - lam / nS) / gammai).flatten() # Otherwise, shrink S following eq.(9)
        
        # compute the gain in the objective funtion at this iteration
        new_norm_beta = np.linalg.norm(newbeta)

        gain[i] = (gammai*(norm_beta[i] + new_norm_beta)/2 + lam)*(norm_beta[i]-new_norm_beta) + \
                    S @ (newbeta - beta[i,:])

        if verbose > 0:
            logging.info(f"[BCD] Iter {itc}, update block {asi}, gain={gain[i]}\n")

        # update beta
        beta[i, :] = newbeta.copy()
        norm_beta[i] = new_norm_beta

        itc += 1
            
    return beta
    
    
def multiplyXtXbysparse(n, active_set, beta, w):
    """"""
    a, p = beta.shape
    
    if a == 0:
        return np.zeros((n-1, p))
    else:
        # First multiply beta by the weights
        beta = np.multiply(beta, np.array([w[active_set]]*p).transpose().reshape(-1, p))
        # compute the matrix s of increments of r
        s = np.zeros((n-1, p))
        s[active_set.flatten(),:] = beta.copy()
        s = np.cumsum(s[::-1,:], 0)[::-1,:]
        u = (active_set.flatten() + 1) @ beta
        s = s - np.array([u]*(n-1)) / n
        r = np.cumsum(s, 0)
        # then multiply the rows by the weights
        r = np.multiply(r, np.array([w]*p).transpose())

    return r


def reconstruct_signal(n, active_set, beta, weights, meansignal=None):
    """ Reconstruct a piecewise-constant profile from its increments."""
    
    a, p = beta.shape
    
    if meansignal is None:
        meansignal = np.zeros((1,p))
    
    signal = np.zeros((n, p))
    signal[active_set.flatten()+1, :] = np.multiply(beta, np.array([weights[active_set.flatten()]]*p).transpose())
    signal = np.cumsum(signal, 0)
    
    # Set the mean signal to meansignal
    m = meansignal - np.mean(signal, 0)
    signal = signal + np.array([m]*n)
    
    return signal


def optimize(beta, active_set, lam, xty, y, w, verbose=1):
    """Solve the group fused Lasso optimization problem"""

    tol = 1e-8
    n, p = xty.shape
    n += 1
    solved = False
    
    t_start = time.time()
    
    itc = 1
    meansignal = np.mean(y, axis=0)
    res = {}
    res['obj'] = np.zeros(10000)
    
    while not solved:
        # optimize gfl given current active set
        beta = block_coordinate_descent(beta, active_set, xty[active_set, :], n, w, lam, verbose=verbose-1)

        # update active set
        
        if len(beta) == 0:
            S = -xty # dim: n-1 x p
            normS = np.sum(S**2, axis=1) # dim n-1
        else:
            nonzero_coef = (np.sum(beta**2, axis=1) !=0)
            active_set = active_set[nonzero_coef,:]
            beta = beta[nonzero_coef, :]

            # Check optimality
            S = multiplyXtXbysparse(n, active_set, beta, w) - xty
            normS = np.sum(S**2, axis=1)
        
        if len(active_set) == 0:
            maxnormS = np.max(normS)
            imax = np.argmax(normS)
            if maxnormS < lam**2 + tol:
                solved = True
            else:
                # this should only occur in the first iteration
                active_set = np.array([[imax]])
                beta = np.zeros((1, p))
        else:
            # for i in AS: normS[i] = lam^2
            # for i not in AS: normS[i] < lam^2
            
            lagr = max(normS[active_set])
            lagr = min(lagr, lam**2)
            nonas = np.setdiff1d(np.arange(n-1), active_set)
            yy = np.max(normS[nonas])
            i = np.argmax(normS[nonas])
            
            if (len(nonas) == 0 or yy < lagr + tol):
                # optimality conditions are met; we found the global solution
                solved = True
                logging.info("Optimality conditions are met!!!")
            else:
                # Otherwise we add the block that violates most the optimality condition
                active_set = np.vstack((active_set, nonas[i]))
                beta = np.vstack((beta, np.zeros((1,p))))
          
        X = reconstruct_signal(n, active_set, beta, w, meansignal)    
        XE = X[1:, :] - X[:-1, :]
        
        res['obj'][itc] = 0.5 * np.linalg.norm(X-y)**2 + lam * onetwo_norm(XE.T)
        res['beta'] = beta
        res['active_set'] = active_set
        res['X'] = X
        
        if verbose >= 1:
            logging.info(f"[optimize] Iteration {itc}\n beta:\n {beta} \n active_set:\n {active_set}\n Obj: {res['obj'][itc]} \n"
            f"-----------------------------------")
        itc += 1
        
    t_end = time.time() - t_start
    logging.info(f"Total time: {t_end}")
    return res
