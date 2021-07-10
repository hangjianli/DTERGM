"""
Generate change statistics from temporal network data -- test
"""
import numpy as np
import argparse
import utils
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula
from rpy2.robjects import numpy2ri
numpy2ri.activate()

network = importr(
    'network',
    on_conflict="warn",
    robject_translations={
        'as.tibble.network': 'as_tibble_network',
        'as_tibble.network': 'as_tibble_network_'
    })

ergm = importr('ergm')
readRDS = robjects.r['readRDS']



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate synthetic data as drawn from an MRF via Gibbs sampling.')

    # I/O settings
    parser.add_argument('outdir', help='The directory to which all of the output will be saved.')
    parser.add_argument('--data_name', help='Name of the raw input data.')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Print detailed progress information to the console.'
                             ' 0=none, 1=high-level only, 2=all details.')

    parser.add_argument('-f', '--form_terms', nargs="+", default=['edges'], help='Formation network statistics.')
    parser.add_argument('-d', '--diss_terms', nargs="+", default=['edges'], help='Dissolution network statistics.')
    # parser.add_argument('--data_format', choices=['real', 'synthetic'], default='synthetic',
    #                     help='The format of input data. Real data is a list of adj matrices while synthetic'
    #                          'data is a list with nw, terms, and extra elements.')
    parser.set_defaults()
    # Get commands from command line
    args = parser.parse_args()
    np.seterr(all='raise')

    # Get the input and output filenames
    output_dir = args.outdir + ('' if args.outdir.endswith('/') else '/')
    args_dir = utils.make_directory(output_dir, 'args')
    print(vars(args))
    utils.save_args(args, args_dir + 'args.txt')
    # input_data = args.data_name + ('' if args.data_name.endswith('.rds') else '.rds')
    # H_outfile = output_dir + args.data_name + "_H_file.txt"
    # y_outfile = output_dir + args.data_name + "_y.txt"
    #
    # # load the data
    # data = readRDS("../data/" + input_data)
    # data = np.array(data).astype(int)
    # t = len(data)
    # n = len(data[0])
    # p = len(args.form_terms) + len(args.diss_terms)
    #
    # print(f"time series length: {t}")
    # print(f"network size: {n} x {n}")
    # print(f"statistics dimension: {p}")
    #
    # H = np.zeros((t, n ** 2, p))  # compute H: T x E x p
    # for i in range(t):
    #     if i == 0:
    #         lr_form, lr_diss = compute_change_statistics(
    #             yt0=data[i],
    #             yt1=data[i + 1],
    #             form_terms=args.form_terms,
    #             diss_terms=args.diss_terms)
    #     else:
    #         lr_form, lr_diss = compute_change_statistics(
    #             yt0=data[i - 1],
    #             yt1=data[i],
    #             form_terms=args.form_terms,
    #             diss_terms=args.diss_terms)
    #     H[i, :, :] = np.concatenate((lr_form, lr_diss), axis=2).reshape(n ** 2, p)  # row-first fill
    #
    # y = data.reshape(t, -1)
    #
    # print('Saving change statistics...')
    # print(f"shape of H: {H.shape}")
    # print(f"shape of y: {y.shape}")
    #
    # np.savetxt(H_outfile, H.reshape(H.shape[0], -1))
    # np.savetxt(y_outfile, y)
    # print("Finished!")