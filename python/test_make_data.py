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
import matplotlib.pyplot as plt
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

def plot_theta(theta, filename, thr=3, change_pts=None):
    """Plot """
    t, p = theta.shape
    fig, axx = plt.subplots(1, 1, figsize=(21, 5))
    norm_diff = np.linalg.norm(np.diff(theta, axis=0), ord=2, axis=1)
    estimated_change_pts = np.arange(t-1)[norm_diff > thr] + 1
    print(f"The estimated change points are at {estimated_change_pts}")
    axx.plot(np.arange(1, t), norm_diff)
    if change_pts is not None:
        for cp in change_pts:
            axx.vlines(x = cp, ymin=0, ymax=max(norm_diff), ls = '--', color = 'r', label='True change points')
    # for cp in estimated_change_pts:
    #     axx.vlines(x = cp, ymin=0, ymax=max(norm_diff), ls = '--', color = 'g', label='Est. change points')
    axx.set_title("l2-norm difference in theta between t and t+1")
    fig.tight_layout()
    plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    simulate_random_graphs(10)