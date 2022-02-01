import numpy as np
import utils
from gfl import *
from scipy.linalg import block_diag
from functools import partial
import time


class STERGMGraph:
    def __init__(
        self,
        lam,
        admm_alpha=100,
        rel_tol=1e-7,
        max_steps=100,
        newton_max_steps=100,
        converge_tol=1e-7,
        gd_lr=0.001,
        gd_epochs=500,
        gfl_tol=1e-6,
        gfl_maxit=1000,
        solver="newton",
        verbose=1
    ):

        self.admm_alpha = admm_alpha
        self.rel_tol = rel_tol
        self.converge_tol = converge_tol
        self.max_steps = max_steps
        self.newton_max_steps = newton_max_steps
        self.verbose = verbose
        self.lam = lam
        self.lr = gd_lr
        self.gfl_tol = gfl_tol,
        self.gfl_maxit = gfl_maxit
        self.epochs = gd_epochs
        self.solver = solver
        self.mu = None

    def load_data(self, h, y, t, p):
        self.H = h #  T x E x p
        self.y = y
        self.t = t
        self.p = p


    def mple_solution_path(self, lam_range=[0.01, 100], bins=50, solver='newton'):
        lambdas = np.exp(np.linspace(np.log(lam_range[1]), np.log(lam_range[0]), bins))
        bic = np.zeros(bins)
        best_lambda = None
        theta_warmstart = None
        for i, lam in enumerate(lambdas):
            self.lam = lam
            pre_res = theta_warmstart
            print(f"[INFO] lambda = {lam}")


        res = self.mple(initial_values=pre_res, solver=solver)
        # self.compute_bic(res)
        # res['bic']
        # to be continued

    def mple(self, initial_values=None, solver='newton', tau_inc=2, tau_dec=2, m=10):
        """Maximum pseudo-likelihood estimation of theta via ADMM"""

        if initial_values:
            theta = initial_values['theta']
            z = initial_values['z']
            u = initial_values['u']
        else:
            theta = np.random.normal(loc=0, scale=0.2, size=(self.t, self.p))
            u = np.zeros((self.t, self.p))
            z = np.random.normal(loc=0, scale=0.2, size=(self.t, self.p))

        converged = self.converge_tol + 1.0
        steps = 0

        # while not converged
        while converged > self.converge_tol and steps < self.max_steps:
            if self.verbose:
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

            if self.verbose:
                print("[INFO] Updating z...")

            # update z by solving group fused lasso (MATLAB)
            z_old = np.copy(z)
            U = theta + u  # observed singal, T x p

            Y = matlab.double(U.tolist())
            _lam = matlab.double([self.lam])
            _tol = matlab.double([self.gfl_tol])
            _maxit = matlab.int16([self.gfl_maxit])

            start = time.time()
#             res = eng.admm_gfl(Y, _lam, _tol, _maxit)
            res = gflasso(Y, _lam)
            end = time.time()
            if self.verbose:
                print(f"[INFO] Group fused Lasso converged in: {end - start}.")
            z = np.asarray(res['_lam']['X']).squeeze().T  # T x p # TODO
            assert z.shape == (self.t, self.p)

            dual_residual = z - z_old
            primal_residual = theta - z

            if self.verbose > 0:
                print("[INFO] Updating u...")

            # update u
            u += primal_residual

            primal_resnorm = np.sqrt(np.mean(primal_residual.flatten() ** 2))
            dual_resnorm = np.sqrt(np.mean(dual_residual.flatten() ** 2))

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

        return {'theta': theta, 'z': z, 'u': u}


    def compute_bic(self, est_chapts):
        pass






    def theta_update(self, z, u, theta):
        """Update theta via Newton method"""

        g = partial(self.theta_update_grad_f, z, u)
        h = partial(self.theta_update_hess_f)

        converged = self.rel_tol + 1.
        steps = 0

        while converged > self.rel_tol and steps < self.newton_max_steps:
            if self.verbose > 1:
                print(f"[INFO] Inner Step # {steps}. Current diff is {converged}")

            hess_f = h(theta)
            grad_f = g(theta)

            # delta_theta = (-np.linalg.inv(hess_f) @ grad_f).reshape(theta.shape)
            # diagonal hessian
            hess = np.diagonal(hess_f)
            delta_theta = (-grad_f / hess).reshape(theta.shape)

            converged = np.abs(np.sum(delta_theta.flatten()) / 2.)
            if converged <= self.rel_tol:
                break

            theta += delta_theta

            if self.verbose > 1:
                print(f"[INFO] Loss = {self.loss_lr(theta, z, u)}")
                print(f"max of mu is {np.max(self.mu)}")
            steps += 1

        return theta  # T x p

    def theta_update_hess_f(self, theta):
        """Update the Hessian in Newton's step

        theta: T x p

        """
        Z = np.sum(self.H * theta[:, np.newaxis, :], axis=2).squeeze()  # compute log odds Zhat: T x E x 1
        if self.verbose > 2:
            print(f"[INFO] Z(TxE): \n {utils.pretty_str(Z, 3)}")
            print(f"[INFO] H(TExTp): \n {utils.pretty_str(self.H, 3)}")
        mu = sigmoid(Z)  # mu: T x E
        W = mu * (1 - mu)  # W: T x E
        _hess1 = self.H.transpose(0, 2, 1) * W[:, np.newaxis, :]
        _hess2 = np.zeros((self.t, self.p, self.p))
        for i in range(self.t):
            _hess2[i, :, :] = _hess1[i, :, :] @ self.H[i, :, :]
        hess = block_diag(*_hess2) + self.admm_alpha * np.identity(self.p * self.t)
        self.mu = mu

        return hess  # Tp x Tp

    def theta_gd_update(self, z, u, theta):
        """Update theta via gradient descent """
        loss = []
        g = partial(self.theta_update_grad_f, z, u)

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

    def theta_update_grad_f(self, z, u, theta):
        """Update the gradient in Newton step """

        Z = np.sum(self.H * theta[:, np.newaxis, :], axis=2).squeeze()  # compute log odds Zhat: T x E x 1
        self.mu = sigmoid(Z)  # mu: T x E
        if self.verbose > 1:
            print(f"[INFO] max mu : \n {np.max(self.mu)}")
        # y = self.X.reshape(self.t, -1)  # T x E
        pnlty = self.admm_alpha * (theta - z + u)  # T x p
        grad = - np.sum(self.H * (self.y - self.mu)[:, :, np.newaxis], axis=1).squeeze() + pnlty

        return grad.flatten()  # 1d Tp

    def loss_lr(self, theta, z, u):
        """First objective function (LR) in ADMM"""
        # y = self.X.reshape(self.t, -1)  # T x E
        # negative pseudo log-likelihood
        predict_1 = self.y * np.log(self.mu)
        predict_0 = (1 - self.y) * np.log(1 - self.mu)
        penalty = np.sum((theta - z + u) ** 2)

        return -np.sum(predict_1 + predict_0) + self.admm_alpha / 2 * penalty


def compute_mle():
    pass


def sigmoid(x):
    return 1 / (1 + safe_exp(-x))


def safe_exp(x):
    return np.exp(x.clip(-50., 50.))
