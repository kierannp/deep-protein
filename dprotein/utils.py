from glob import glob
import numpy as np
import pandas as pd
import subprocess
import os

def get_mean_descriptor(path, remove_time = 0 , descriptor = 'rmsd'):
    try:
        if descriptor == 'rg':
            fname = '/gyrate.xvg'
        if descriptor == 'rmsd':
            if not os.path.exists(path + '/rmsd.xvg'):
                os.chdir(path)
                subprocess.run('module load gromacs && echo echo "1 1 " | gmx rms -s npt.gro -f production.xtc -o rmsd.xvg', shell=True)
                os.chdir('-')
            fname = '/rmsd.xvg'
        else:
            raise Exception('Descriptor not implemented!')
        
        with open(path + fname) as f:
            lines = f.readlines()
            f.close()
        values = []
        start = 0
        for i, l in enumerate(lines):
            if descriptor == 'rg':
                if "s3 legend" in l:
                    start = i
                    break
            if descriptor == 'rmsd':
                if "@ subtitle " in l:
                    start = i
                    break
        for l in lines[start+1:]:
            temp = list(map(float, l.split()[1:]))
            values.append(temp[0])
        return np.mean(np.array(values)[remove_time*1000:])
    except:
        return None

def process_results(sim_path, descriptor = 'rmsd', remove_time = 0):
    paths = glob(sim_path+'/*')
    values, donors, acceptors, h2o = [], [], [], []
    for p in paths:
        (donor_percent, h2o_percent, _) = p.split('/')[-1].split('_')
        donor_percent, h2o_percent = float(donor_percent), float(h2o_percent)
        des_percent = 1.0-h2o_percent
        acceptor_percent = (1.0 - donor_percent)*des_percent
        donor_percent = donor_percent*des_percent

        donors.append(donor_percent)
        acceptors.append(acceptor_percent)
        h2o.append(h2o_percent)

        values.append(get_mean_descriptor(p, descriptor = descriptor, remove_time=remove_time))
    return pd.DataFrame({'hbd':donors, 'hba':acceptors, 'h2o':h2o, descriptor:values})