from glob import glob
import numpy as np
import pandas as pd

def get_mean_gyration(path):
    try:
        with open(path + '/gyrate.xvg') as f:
            lines = f.readlines()
            f.close()
        gyrations = []
        start = 0
        for i, l in enumerate(lines):
            if "s3 legend" in l:
                start = i
                break
        for l in lines[start+1:]:
            values = list(map(float, l.split()[1:]))
            gyrations.append(values[1])
        
        return np.mean(np.array(gyrations))
    except:
        return None
def process_results(sim_path):
    paths = glob(sim_path+'/*')
    gyrations, donors, acceptors, h2o = [], [], [], []
    for p in paths:
        (donor_percent, h2o_percent, rep) = p.split('/')[-1].split('_')
        donor_percent, h2o_percent = float(donor_percent), float(h2o_percent)
        des_percent = 1.0-h2o_percent
        acceptor_percent = (1.0 - donor_percent)*des_percent
        donor_percent = donor_percent*des_percent

        donors.append(donor_percent)
        acceptors.append(acceptor_percent)
        h2o.append(h2o_percent)
        gyrations.append(get_mean_gyration(p))
    return pd.DataFrame({'hbd':donors, 'hba':acceptors, 'h2o':h2o, 'rg':gyrations})