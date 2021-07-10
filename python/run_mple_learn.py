import numpy as np
import mple_learn_new
import matlab.engine
import argparse

eng = matlab.engine.start_matlab()
eng.addpath(eng.genpath('../gfl/'), nargout=0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run MPLE algorithm.')

    # I/O settings
    parser.add_argument('outdir', help='The directory to which all of the output will be saved.')


    # mple settings
    parser.add_argument('--lam', default=10, help="")
    parser.add_argument('--admm_alpha', default=1, help="Initial ADMM step size.")
    parser.add_argument('--solver', type=str, default='newton', help="Choose the solver for logistic regression.")
    parser.add_argument('--rel_tol', type=float, default=1e-7, help="Stopping criterion for change in theta.")
    parser.add_argument('--max_steps', type=int, default=100, help="Max iteration for ADMM.")
    parser.add_argument('--max_steps_newton', type=int, default=100, help="Max iterations for the Newton step.")
    parser.add_argument('--conv_tol', type=float, default=1e-5, help="Stopping criterion for max of dual/primal residuals.")
    parser.add_argument('--gd_lr', type=float, default=0.01, help="Gradient descent learning rate")
    parser.add_argument('--gd_epochs', type=int, default=100, help="Number of epochs for gradient descent.")
    parser.add_argument('--gfl_tol', type=float, default=1e-5, help="Stop criterion for group fused lasso.")
    parser.add_argument('--gfl_maxit', type=int, default=1000, help="Max iter for group fused lasso.")
    # collect the command line inputs
    parser.set_defaults()
    args = parser.parse_args()

    np.seterr(all='raise')

    # Load setup info and data
    with open(args.outdir + '/args/args.txt', 'r') as f:
        sim_args = f.read()
    p = len(sim_args['form_terms']) + len(sim_args['diss_terms'])
    H = np.loadtxt(args.outdir + sim_args['data_name'] + '.txt')
    y = np.loadtxt(args.outdir + sim_args['data_name'] + '.txt')
    t = H.shape[0]
    H = H.reshape((t, -1, p))
    n = H.shape[1]
    print(f"(n, p, t) = ({n}, {p}, {t})")





    model = mple_learn_new.STERGMGraph(
        lam=args.lam,
        admm_alpha=args.admm_alpha,
        rel_tol= args.rel_tol,
        max_steps=args.max_steps,
        newton_max_steps=args.max_steps_newton,
        converge_tol=args.conv_tol,
        gd_lr=args.gd_lr,
        gd_epochs=args.gd_epochs,
        gfl_tol=args.gfl_tol,
        gfl_maxit=args.maxit,
        verbose=sim_args["verbose"])

    model.load_data(H, y, t, n, p)

    res = model.mple()
    theta = res['theta']
    u = res['u']
    z = res['z']
    np.savetxt(theta)

    print("Done!")
