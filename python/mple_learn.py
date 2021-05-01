import numpy as np
from scipy.linalg import block_diag

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
network = importr('network', on_conflict="warn",
                  robject_translations={
                      'as.tibble.network': 'as_tibble_network',
                      'as_tibble.network': 'as_tibble_network'
                  }
         )
base = importr('base')
ergm = importr('ergm')
from rpy2.robjects import Formula

class stergmGraph:
    def __init__(
        self,
        X,
        form_terms,
        diss_terms,
        admm_alpha=100,
        rel_tol=1e-6,
        newton_max_steps=100,
        converge_tol=1e-4,
        verbose=True
    ):
        self.X = X
        self.n = len(X[0])
        self.tslen = len(X)
        self.pf = len(form_terms)
        self.pd = len(diss_terms)
        self.admm_alpha = admm_alpha
        self.rel_tol = rel_tol
        self.verbose = verbose
        self.fml_form = Formula('nwf ~ ' + ' + '.join(form_terms))
        self.fml_diss = Formula('nwd ~ ' + ' + '.join(diss_terms))
    
    
    def theta_update_hess(self, theta, X):
        """Update the Hessian in Newton's step 

        theta: T x p
        
        """
        p = self.pf + self.pd
        H = np.zeros((self.tslen, self.n**2, p))  # compute H: T x E x p
        for i in range(self.tslen):
            if i == 0:
                lr_form, lr_diss = self.g_delta(yt0=X[i], yt1=X[i])
            else:
                lr_form, lr_diss = self.g_delta(yt0=X[i-1], yt1=X[i])
            H[i,:,:] = np.concatenate((lr_form, lr_diss), axis=2).reshape(self.n**2, p)   # row-first fill

        Z = np.sum(H * theta[:,np.newaxis, :], axis=2).squeeze() # compute log odds Zhat: T x E x 1 
        mu = sigmoid(Z)         # mu: T x E
        W = mu * (1 - mu)        # W: T x E 
        _hess1 = H.transpose(0,2,1) * W[:,np.newaxis,:] 
        _hess2 = np.zeros((self.tslen, p, p))
        for i in range(self.tslen):
            _hess2[i,:,:] = _hess1[i,:,:] @ H[i,:,:]
        hess = block_diag(*_hess2) + self.admm_alpha * np.identity(p * self.tslen)
        self.H = H
        self.mu = mu
        return hess # Tp x Tp
                
    
    
    def theta_update_grad(self, theta, X, z, u):
        y = X.reshape(self.tslen, -1) # T x E
        pnlty = self.admm_alpha * (theta - z + u) # T x p
        grad = (- self.H.transpose(0,2,1) @ (y - self.mu)).squeeze() + pnlty
        return grad
        
        
    def theta_update(self, theta, X, z, u):
        pass
    
    def g_delta(self, yt0, yt1):
        """ Compute the change stat in mple for all edges  """
        # 1. take union. 2. run ergmMPLE, 3. correct
        ypos = np.logical_or(yt0, yt1).astype(int)
        yneg = np.logical_and(yt0, yt1).astype(int)
        nwpos = network.network(ypos)
        nwneg = network.network(yneg)
        self.fml_form.environment['nwf'] = nwpos
        self.fml_diss.environment['nwd'] = nwneg
        
        # compute change statistic with correction
        lr_form = ergm.ergmMPLE(self.fml_form, output='array')
        lr_form_delta = np.array(lr_form.rx('predictor')).squeeze(axis=0)
        lr_diss = ergm.ergmMPLE(self.fml_diss, output='array')
        lr_diss_delta = np.array(lr_diss.rx('predictor')).squeeze(axis=0)

        # correct and set diagonal to 0
        for i in range(self.pf):
            lr_form_delta[:,:,i][yt0==1] = 0
            np.fill_diagonal(lr_form_delta[:,:,i], 0)
        for i in range(self.pd):
            lr_diss_delta[:,:,i][yt0==0] = 0
            np.fill_diagonal(lr_diss_delta[:,:,i], 0)
        
        return lr_form_delta, lr_diss_delta # n x n x p1 (p2)
        

def sigmoid(x):
    return 1/(1 + np.exp(-x))            





#     def gstat_delta(self, yt0, yt1):
#         """ Compute the change stat in mple for all edges  """
#         gdelta = np.zeros((self.p, self.n**2))
#         for j in range(self.n):
#             for i in range(self.n):
#                 if i != j:
#                     cache = yt1[i,j]
#                     yt1[i,j] = 1
#                     gstat1 = self.gstat(yt0, yt1)
#                     yt1[i,j] = 0
#                     gstat2 = self.gstat(yt0, yt1)
#                     yt1[i,j] = cache
#                     gdelta[:, j*self.n+i] = gstat1 - gstat2
    
#         return gdelta

#     def gstat(self, yt0, yt1):
#         '''Compute the suff stat in stergm at time t'''
#         ypos = np.logical_or(yt0, yt1).astype(int)
#         yneg = np.logical_and(yt0, yt1).astype(int)
        
#         nwpos = network.network(ypos)
#         nwneg = network.network(yneg)
#         self.fml_form.environment['nwf'] = nwpos
#         self.fml_diss.environment['nwd'] = nwneg

#         gpos = base.summary(self.fml_form)
#         gneg = base.summary(self.fml_diss)
        
#         return np.append(gpos, gneg)
    
