"""
p =
H =
y =
"""

import numpy as np
import mple_learn
import matlab.engine
import argparse
import ast
import time
import utils
import matplotlib.pyplot as plt

eng = matlab.engine.start_matlab()
eng.addpath(eng.genpath('../gfl/'), nargout=0)

def plot_theta(theta, filename, thr=3, change_pts=None):
    """Plot """
    t, p = theta.shape
    fig, axx = plt.subplots(1, 1, figsize=(21, 5))
    norm_diff = np.linalg.norm(np.diff(theta, axis=0), ord=2, axis=1)
    estimated_change_pts = np.arange(t-1)[norm_diff > thr] + 1
    print(f"The estimated change points are {estimated_change_pts}")
    axx.plot(np.arange(1, t), norm_diff)
    if change_pts is not None:
        for i, cp in enumerate(change_pts):
            axx.vlines(x = cp, ymin=0, ymax=max(norm_diff), ls = '--', color = 'r',
                       label='True change points' if i == 0 else "")
    axx.set_title("l2-norm difference in theta between t and t+1", fontsize=30)
    fig.tight_layout()
    plt.legend(fontsize=20)
    plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run MPLE algorithm.')

    # I/O settings
    parser.add_argument('outdir', help='The directory to which all of the output will be saved.')
    # mple settings
    parser.add_argument('--lam', type=float, nargs="+", default=[10], help="")
    parser.add_argument('--admm_alpha', type=float, default=1, help="Initial ADMM step size.")
    parser.add_argument('--solver', type=str, default='newton', help="Choose the solver for logistic regression.")
    parser.add_argument('--rel_tol', type=float, default=1e-7, help="Stopping criterion for change in theta.")
    parser.add_argument('--max_steps', type=int, default=100, help="Max iteration for ADMM.")
    parser.add_argument('--max_steps_newton', type=int, default=100, help="Max iterations for the Newton step.")
    parser.add_argument('--conv_tol', type=float, default=1e-5,
                        help="Stopping criterion for max of dual/primal residuals.")
    parser.add_argument('--gd_lr', type=float, default=0.01, help="Gradient descent learning rate")
    parser.add_argument('--gd_epochs', type=int, default=100, help="Number of epochs for gradient descent.")
    parser.add_argument('--gfl_tol', type=float, default=1e-5, help="Stop criterion for group fused lasso.")
    parser.add_argument('--gfl_maxit', type=int, default=1000, help="Max iter for group fused lasso.")

    parser.add_argument('--true_cp', type=int, nargs="+", default=None, help="True change points (for plot).")
    parser.add_argument('--cp_thr', type=float, default=3, help="Threshold for estimating change points.")
    # collect the command line inputs

    parser.set_defaults()
    args = parser.parse_args()

    np.seterr(all='raise')

    # Load setup info and data
    print("Loading network statistics and edge data...")
    experiment_dir = args.outdir + ('' if args.outdir.endswith('/') else '/')
    with open(experiment_dir + 'args/args.txt', 'r') as f:
        sim_args_str = f.read()
    sim_args = ast.literal_eval(sim_args_str)
    p = len(sim_args['form_terms']) + len(sim_args['diss_terms'])
    H = np.loadtxt(experiment_dir + sim_args['data_name'] + '_H.txt')
    y = np.loadtxt(experiment_dir + sim_args['data_name'] + '_y.txt')
    t = H.shape[0]
    H = H.reshape((t, -1, p)) # t x n^2(E) x p
    n = np.sqrt(H.shape[1]).astype(int)
    print(f"Data has dimension (t, n, p): ({t}, {n}, {p})")

    # Get the output filenames
    result_dir = utils.make_directory(experiment_dir, 'results')
    args_dir = utils.make_directory(experiment_dir, 'args')
    utils.save_args(args, args_dir + 'args_model.txt')

    theta_outfile = result_dir + 'theta_' + sim_args['data_name'] + ".txt"
    u_outfile = result_dir + 'u_' + sim_args['data_name'] + ".txt"
    z_outfile = result_dir + 'z_' + sim_args['data_name'] + ".txt"
    theta_plot_dir = result_dir + 'est_theta_diff.png'

    print('Initialize STERGM model...')
    model = mple_learn.STERGMGraph(
        lam=args.lam,
        admm_alpha=args.admm_alpha,
        rel_tol=args.rel_tol,
        max_steps=args.max_steps,
        newton_max_steps=args.max_steps_newton,
        converge_tol=args.conv_tol,
        gd_lr=args.gd_lr,
        gd_epochs=args.gd_epochs,
        gfl_tol=args.gfl_tol,
        gfl_maxit=args.gfl_maxit,
        verbose=sim_args["verbose"])

    model.load_data(H, y, t, p)

    print("Run mple to estimate theta...")
    start = time.time()
    res = model.mple(solver=args.solver)
    end = time.time()
    print(f"[INFO] MPLE finished in: {end - start}.")
    theta = res['theta']
    u = res['u']
    z = res['z']

    print("Saving estimated parameters...")
    np.savetxt(theta_outfile, theta)
    np.savetxt(u_outfile, u)
    np.savetxt(z_outfile, z)

    plot_theta(theta, theta_plot_dir, args.cp_thr, args.true_cp)
    print("Done!")
