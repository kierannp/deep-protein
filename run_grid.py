import sys
sys.path.insert(1, '~/projects/deep-protein')
from dprotein import *
from dprotein.utils import *
import argparse



parser = argparse.ArgumentParser(
                    prog='ternary_grid_search',
                    description='runs grid search over MD simulation of a ternary phase diagram',
                    )

parser.add_argument('--n_solvent_mols', type=int, nargs='+', default=2250)
parser.add_argument('-c','--mol_charges', nargs='+', default= [1, -1, 0, -4])
parser.add_argument('--worker_num', type=int)  
parser.add_argument('--workers', nargs=1, type=int, default=4)
# parser.add_argument('-c', '--count')      # option that takes a value

args = parser.parse_args()
print(args)


repo_path = '/raid6/homes/kierannp/projects/deep-protein/'

gyration_data = {}

simulation = Deep_Eutectic_Search(
    pdb_file = repo_path + 'data/pdbs/proteins/1OIL_A_no_CAL.pdb',
    hba_file_1 = repo_path + 'data/pdbs/HBA/choline.pdb',
    hba_file_2 = repo_path + 'data/pdbs/HBA/cla.pdb',
    hbd_file = repo_path + 'data/pdbs/HBD/urea.pdb',
    h2o_file = repo_path + 'data/pdbs/water.pdb',
    working_dir = repo_path + 'sandbox',
    charges = {'hba_1' : args.mol_charges[0], 
               'hba_2' : args.mol_charges[1], 
               'hbd' : args.mol_charges[2], 
               'protein' : args.mol_charges[3]
               },
    total_solvent_mols = args.n_solvent_mols
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
            buffer_size=3.0,
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

# create grid on ternary diagram
n = 13
nx, ny = (n, n)
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
xv, yv = np.meshgrid(x, y)
z = np.ones(n**2) - xv.flatten() - yv.flatten()

data = {
    'hbd':xv.flatten()[z>=0],
    'hba':z[z>=0],
    'h2o':yv.flatten()[z>=0],
}

np.seterr(invalid='ignore') # needed for divide by zero in np.where

test_points = {
    'percent_hbd' : np.where(data['hbd']+data['hba'] == 0, 0, data['hbd'] / (data['hbd']+data['hba'])), 
    'percent_h2o' : data['h2o']
}
test_points = [list(p) for p in zip(test_points['percent_hbd'], test_points['percent_h2o'])]


# get previous simulations
fnames = glob('./sandbox/*')
hba, h2o = [], []
for f in fnames:
    file_info = f.split('/')[-1].split('_')
    hba.append(float(file_info[0]))
    h2o.append(float(file_info[1]))

# remove sampled simulations from grid
test_points = np.array(test_points)
sampled_points = []
for i, (a, b) in enumerate(zip(hba,h2o)):
    for j in range(test_points.shape[0]):
        if np.allclose(np.array([a,b]), test_points[j,:]):
            sampled_points.append(j)

# split grid into sections for parallel search
clean_grid = np.delete(test_points, np.array(list(set(sampled_points))), axis=0)
n_points = clean_grid.shape[0]
if args.worker_num <= 1:
    section = clean_grid[:n_points//args.workers]
else:
    section = clean_grid[(args.worker_num - 1) * (n_points//args.workers) : args.worker_num * (n_points//args.workers)]


for point in section:
    run_simulation(point)
