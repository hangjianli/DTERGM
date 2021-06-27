import numpy as np
from scipy.linalg import block_diag
from functools import partial
from utils import *
import matlab.engine
import time

eng = matlab.engine.start_matlab()
eng.addpath(eng.genpath('../gfl/'), nargout=0)

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
        lam,
        admm_alpha=100,
        rel_tol=1e-7,
        max_steps=100,
        newton_max_steps=100,
        converge_tol=1e-7,
        gd_lr = 0.001,
        gd_epochs = 500,
        gfl_tol = 1e-6,
        gfl_maxit=1000,
        verbose=1
    ):
        self.X = X
        self.n = len(X[0])
        self.tslen = len(X)
        self.pf = len(form_terms)
        self.pd = len(diss_terms)
        self.p = len(form_terms) + len(diss_terms)
        self.admm_alpha = admm_alpha
        self.rel_tol = rel_tol
        self.converge_tol = converge_tol
        self.max_steps = max_steps
        self.newton_max_steps = newton_max_steps
        self.verbose = verbose
        self.fml_form = Formula('nwf ~ ' + ' + '.join(form_terms))
        self.fml_diss = Formula('nwd ~ ' + ' + '.join(diss_terms))
        self.lamseq = lam
        self.lr = gd_lr
        self.gfl_tol = gfl_tol,
        self.gfl_maxit = gfl_maxit
        self.epochs = gd_epochs
        self.mu = None
        self.H = np.zeros((self.tslen, self.n**2, self.p))  # compute H: T x E x p
        
        for i in range(self.tslen):
            if i == 0:
                lr_form, lr_diss = self.g_delta(yt0=self.X[i], yt1=self.X[i+1])
            else:
                lr_form, lr_diss = self.g_delta(yt0=self.X[i-1], yt1=self.X[i])
            self.H[i,:,:] = np.concatenate((lr_form, lr_diss), axis=2).reshape(self.n**2, self.p)   # row-first fill
    
    
    def mple(self, 
             initial_values=None, 
             solver = 'newton', # 'newton' or 'gd'
             tau_inc=2, 
             tau_dec=2, 
             m=10):
        """Maximum pseudo-likelihood estimation of theta via ADMM"""
        
        if initial_values:
            theta = initial_values['theta']
            z = initial_values['z']
            u = initial_values['u']
        else:
            theta = np.random.normal(size=(self.tslen, self.p))
            z = np.zeros(size=(self.tslen, self.p))
            u = np.zeros(size=(self.tslen, self.p))
        
        converged = self.converge_tol + 1.0
        steps = 0
        
        # while not converged
        while converged > self.converge_tol and steps < self.max_steps:
            if self.verbose > 0:
                print(f"[INFO] ADMM step #{steps}")
                print("[INFO] Updating theta...")
            
            # update theta
            if solver == 'newton':
                start = time.time()
                theta = self.theta_update(z, u, theta)
                end = time.time()
                if self.verbose:
                    print(f"[INFO] Newton converged in: {end - start}.")
            else:
                start = time.time()
                theta = self.theta_gd_update(z, u, theta)
                end = time.time()
                if self.verbose:
                    print(f"[INFO] Gradient descent converged in: {end - start}.")
                
            
            if self.verbose > 0:
                print("[INFO] Updating z...")
            
            # update z by solving group fused lasso (MATLAB)
            z_old = np.copy(z)
            U = theta + u # observed singal, T x p
            
            Y = matlab.double(U.tolist())
            _lam = matlab.double([lam * self.admm_alpha for lam in self.lamseq])
            _tol = matlab.double([self.gfl_tol])
            _maxit = matlab.int16([self.gfl_maxit])

            start = time.time()
            res = eng.admm_gfl(Y, _lam, _tol, _maxit)  
            end = time.time()
            if self.verbose:
                    print(f"[INFO] Group fused Lasso converged in: {end - start}.")
            z = np.asarray(res['X']).squeeze().T # T x p
            assert z.shape == (self.tslen, self.p)
            
            dual_residual = z - z_old
            primal_residual = theta - z
            
            if self.verbose > 0:
                print("[INFO] Updating u...")
                
            # update u    
            u += primal_residual
            
            primal_resnorm = np.sqrt(np.mean(primal_residual.flatten()**2))
            dual_resnorm = np.sqrt(np.mean(dual_residual.flatten()**2))
            
            converged = max(primal_resnorm, dual_resnorm)
            
            if primal_resnorm > dual_resnorm * m:
                self.admm_alpha *= tau_inc
                u /= tau_inc
                if self.verbose > 0:
                    print("-------------------------------------------------")
                    print(f"[INFO] increasing alpha to {self.admm_alpha}")
            elif dual_resnorm > primal_resnorm * m:
                self.admm_alpha /= tau_dec
                u *= tau_dec
                if self.verbose > 0:
                    print("-------------------------------------------------")
                    print(f"[INFO] decreasing alpha to {self.admm_alpha}")
                
            if self.verbose > 0:
                print(f"[INFO] max mu :  {np.max(self.mu)}")
                print(f"[INFO] dual_resnorm: {dual_resnorm:.6f}")
                print(f"[INFO] primal_resnorm: {primal_resnorm:.6f}")
                print(f"[INFO] convergence: {converged:.6f}")
                print("-------------------------------------------------")
            
            steps += 1
    
        if self.verbose:
            print("[INFO] ADMM finished!")
    
        return theta, z, u, converged
    
   
    
    def theta_update(self, z, u, theta):
        """Update theta via Newton method"""
        
        g = partial(self.theta_update_grad_f, z, u)
        h = partial(self.theta_update_hess_f)
        
        converged = self.rel_tol + 1.
        steps = 0
        
        while converged > self.rel_tol and steps < self.newton_max_steps:
            if self.verbose > 1:
                print(f"[INFO]\tInner Step # {steps}. Current diff is {converged}")
                
            hess_f = h(theta)            
            grad_f = g(theta)
  
            delta_theta = (-np.linalg.inv(hess_f) @ grad_f).reshape(theta.shape)

            converged = np.abs(np.sum(delta_theta.flatten()) / 2.)
            if converged <= self.rel_tol:
                break
                
            theta += delta_theta
            
            if self.verbose > 1:
                print(f"[INFO] Loss = {self.loss(theta, z, u)}")
                print(f"max of mu is {np.max(self.mu)}")
            steps += 1

        return theta # T x p
    
    
    def theta_update_hess_f(self, theta):
        """Update the Hessian in Newton's step 

        theta: T x p
        
        """
        Z = np.sum(self.H * theta[:,np.newaxis, :], axis=2).squeeze() # compute log odds Zhat: T x E x 1 
        if self.verbose > 2:
            print(f"[INFO] Z(TxE): \n {pretty_str(Z, 3)}")
            print(f"[INFO] H(TExTp): \n {pretty_str(H, 3)}")
        mu = sigmoid(Z)         # mu: T x E
        W = mu * (1 - mu)        # W: T x E 
        _hess1 = self.H.transpose(0,2,1) * W[:,np.newaxis,:] 
        _hess2 = np.zeros((self.tslen, self.p, self.p))
        for i in range(self.tslen):
            _hess2[i,:,:] = _hess1[i,:,:] @ self.H[i,:,:]
        hess = block_diag(*_hess2) + self.admm_alpha * np.identity(self.p * self.tslen)
        self.mu = mu
        
        return hess # Tp x Tp
              
        
    def theta_gd_update(self, z, u, theta):
        """Update theta via gradient descent """
        loss = []
        g = partial(self.theta_update_grad_f, z, u)
        converged = 1.
        
        for it in range(self.epochs):        
            delta_theta = self.lr * g(theta).reshape(theta.shape)
            converged = np.abs(np.sum(delta_theta.flatten()) / 2.)
            if converged < self.rel_tol:
                if self.verbose:
                    print(f"Gradient descent converged. Epoch number {it}")
                break
            theta -= delta_theta
            curloss = self.loss_lr(theta, z, u)
            if self.verbose > 2:
                print(f"[INFO] Loss: {curloss: .5f}")
            loss.append(curloss) 
        self.loss = loss
        
        return theta
        
    def theta_update_grad_f(self,  z, u, theta):
        """Update the gradient in Newton step """

        Z = np.sum(self.H * theta[:,np.newaxis, :], axis=2).squeeze() # compute log odds Zhat: T x E x 1 
        self.mu = sigmoid(Z)         # mu: T x E
        if self.verbose > 1:
            print(f"[INFO] max mu : \n {np.max(self.mu)}")
        y = self.X.reshape(self.tslen, -1) # T x E
        pnlty = self.admm_alpha * (theta - z + u) # T x p
        grad = - np.sum( self.H * (y - self.mu)[:,:,np.newaxis], axis=1).squeeze() + pnlty
        
        return grad.flatten() # 1d Tp
        
    
    def g_delta(self, yt0, yt1):
        """ Compute the change stat in mple for all edges  """
        # 1. take union. 2. run ergmMPLE, 3. adjust
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
        
        
    def loss_lr(self, theta, z, u):
        """First objective function (LR) in ADMM"""
        y = self.X.reshape(self.tslen, -1) # T x E
        # negative pseudo log-likelihood
        predict_1 = y * np.log(self.mu)
        predict_0 = (1 - y) * np.log(1 - self.mu)
        penalty = np.sum((theta - z + u)**2)
        
        return -np.sum(predict_1 + predict_0) + self.admm_alpha / 2 * penalty
        

        
def compute_mle():
    pass
    

def sigmoid(x):
    return 1/(1 + safe_exp(-x))            


def safe_exp(x):
    return np.exp(x.clip(-50.,50.))


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
    