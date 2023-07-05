import sys
# from bayes_opt import BayesianOptimization
# from sklearn.gaussian_process import GaussianProcessRegressor
# from sklearn.gaussian_process.kernels import RBF
# from bayes_opt import BayesianOptimization
# from sklearn.gaussian_process.kernels import Matern, RationalQuadratic, DotProduct, ConstantKernel

from sklearn.base import clone
from skopt import gp_minimize
from skopt.learning import GaussianProcessRegressor
from skopt.learning.gaussian_process.kernels import ConstantKernel, Matern, RationalQuadratic, DotProduct


sys.path.insert(1, '~/projects/deep-protein')
from dprotein import *
from dprotein.utils import *

repo_path = '/raid6/homes/kierannp/projects/deep-protein/'

# Bounded region of parameter space
pbounds = {'donor_percent': (0, 1), 'h2o_percent': (0, 1)}
molar_weights = {'chcl': 139.62, 'urea': 60.06, 'water':18.02} # g/mol
gyration_data = {}

simulation = Deep_Eutectic_Search(
    pdb_file = repo_path + 'data/pdbs/proteins/1OIL_A_no_CAL.pdb',
    hba_file_1 = repo_path + 'data/pdbs/HBA/choline.pdb',
    hba_file_2 = repo_path + 'data/pdbs/HBA/cla.pdb',
    hbd_file = repo_path + 'data/pdbs/HBD/urea.pdb',
    h2o_file = repo_path + 'data/pdbs/water.pdb',
    working_dir = repo_path + 'sandbox',
    charges = {'hba_1' : 1, 'hba_2' : -1, 'hbd' : 0, 'protein' : -4},
    total_solvent_mols = 1750
)

def run_simulation(percents = []):
    n_replicates = 2
    donor_percent, h2o_percent = percents[0], percents[1]
    results = []
    if os.path.exists('/raid6/homes/kierannp/projects/deep-protein/sandbox/{}_{}_0'.format(donor_percent, h2o_percent)):
        return get_mean_gyration('/raid6/homes/kierannp/projects/deep-protein/sandbox/{}_{}_0'.format(donor_percent, h2o_percent))
    for r in range(n_replicates):
        simulation.build_system(
            donor_percent = donor_percent,
            h2o_percent = h2o_percent,
            seed = r
        )
        simulation.apply_forcefield()
        simulation.run_sim(
            em_mdp = repo_path + 'mdps/em.mdp',
            nvt_mdp = repo_path + 'mdps/nvt.mdp',
            npt_mdp = repo_path + 'mdps/npt_373.mdp',
            production_mdp = repo_path + 'mdps/production_5ns.mdp'
        )
        results.append(simulation.analyze_trajectory())
    return np.mean(np.array(results))

df = process_results(sim_path = '/raid6/homes/kierannp/projects/deep-protein/sandbox')
df = df.assign(unnorm_hbd = lambda x: x['hbd']/(1-x['h2o']))
df = df.dropna()


def plot_approximation(gpr, X, Y, X_sample, Y_sample, X_next=None, show_legend=False):
    mu, std = gpr.predict(X, return_std=True)
    plt.fill_between(X.ravel(), 
                     mu.ravel() + 1.96 * std, 
                     mu.ravel() - 1.96 * std, 
                     alpha=0.1) 
    plt.plot(X, Y, 'y--', lw=1, label='Noise-free objective')
    plt.plot(X, mu, 'b-', lw=1, label='Surrogate function')
    plt.plot(X_sample, Y_sample, 'kx', mew=3, label='Noisy samples')
    if X_next:
        plt.axvline(x=X_next, ls='--', c='k', lw=1)
    if show_legend:
        plt.legend()

def plot_acquisition(X, Y, X_next, show_legend=False):
    plt.plot(X, Y, 'r-', lw=1, label='Acquisition function')
    plt.axvline(x=X_next, ls='--', c='k', lw=1, label='Next sampling location')
    if show_legend:
        plt.legend()    
        
def plot_convergence(X_sample, Y_sample, n_init=2):
    plt.figure(figsize=(12, 3))

    x = X_sample[n_init:].ravel()
    y = Y_sample[n_init:].ravel()
    r = range(1, len(x)+1)
    
    x_neighbor_dist = [np.abs(a-b) for a, b in zip(x, x[1:])]
    y_max_watermark = np.maximum.accumulate(y)
    
    plt.subplot(1, 2, 1)
    plt.plot(r[1:], x_neighbor_dist, 'bo-')
    plt.xlabel('Iteration')
    plt.ylabel('Distance')
    plt.title('Distance between consecutive x\'s')

    plt.subplot(1, 2, 2)
    plt.plot(r, y_max_watermark, 'ro-')
    plt.xlabel('Iteration')
    plt.ylabel('Best Y')
    plt.title('Value of best selected sample')

bounds = [(0.0, 1.0), (0.0, 1.0)]
noise  = 0.2
# X_init = np.array([[-0.9], [1.1]])

X_init = df.loc[:,['unnorm_hbd','h2o']].values
Y_init = df.loc[:,['rg']].values

# Use custom kernel and estimator to match previous example
m52 = RationalQuadratic(alpha=1) + DotProduct()
gpr = GaussianProcessRegressor(kernel=m52)

r = gp_minimize(run_simulation,
                bounds,
                base_estimator=gpr,
                # acq_func='LCB',      # Lower confidence bound
                # xi=0.01,            # exploitation-exploration trade-off
                n_calls=10,         # number of iterations
                n_random_starts=1,  # initial samples are provided
                # n_initial_points = - X_init.shape[0],
                # x0=X_init.tolist(), # initial samples
                # y0=Y_init.ravel().tolist(),
)

# Fit GP model to samples for plotting results
gpr.fit(r.x_iters, -r.func_vals)

# Plot the fitted model and the noisy samples
# plot_approximation(gpr, X, Y, r.x_iters, -r.func_vals, show_legend=True)
