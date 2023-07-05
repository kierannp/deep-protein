import shutil
import os
import sys
import pathlib
import subprocess
import numpy as np

class Deep_Eutectic_Search:
    def __init__(self,  pdb_file, hba_file_1, hba_file_2, hbd_file, h2o_file, working_dir, total_solvent_mols, charges):
        self.pdb_file = pdb_file
        self.hba_file_1 = hba_file_1
        self.hba_file_2 = hba_file_2
        self.hbd_file = hbd_file
        self.h2o_file = h2o_file

        self.charges = charges
        self.total_solvent_mols = total_solvent_mols
        self.home = os.path.expanduser("~")
        self.working_dir = working_dir

    def build_system(self, donor_percent, h2o_percent, buffer_size, seed):
        """
        A function that builds the atomistic system

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.data_dir = self.working_dir + '/{}_{}_{}'.format(donor_percent, h2o_percent, seed)
        os.mkdir(self.data_dir)
        os.chdir(self.data_dir)
        local_pdb_file = pathlib.Path(self.pdb_file).stem

        shutil.copy(self.pdb_file, '{}.pdb'.format(local_pdb_file))

        self.n_h2o = int(self.total_solvent_mols * h2o_percent)
        n_des_mol = self.total_solvent_mols - self.n_h2o
        self.n_hba = int((1 - donor_percent) * n_des_mol)
        self.n_hbd = int(donor_percent * n_des_mol)

        system_charge = self.charges['hbd']*self.n_hbd + self.charges['hba_1']*self.n_hba + self.charges['hba_2']*self.n_hba + self.charges['protein']
        if system_charge != 0:
            if system_charge > 0 and self.charges['hba_1'] < 0:
                n_hba_1 = self.n_hba + abs(system_charge // self.charges['hba_1'])
                n_hba_2 = self.n_hba
            elif system_charge > 0 and self.charges['hba_2'] < 0:
                n_hba_1 = self.n_hba
                n_hba_2 = self.n_hba + abs(system_charge // self.charges['hba_2'])
            elif system_charge < 0 and self.charges['hba_1'] > 0:
                n_hba_1 = self.n_hba + abs(system_charge // self.charges['hba_1'])
                n_hba_2 = self.n_hba
            elif system_charge < 0 and self.charges['hba_2'] > 0:
                n_hba_1 = self.n_hba
                n_hba_2 = self.n_hba + abs(system_charge // self.charges['hba_2'])
            else:
                raise Exception('Net system charge, and cannot neutralize')
        else:
            n_hba_1 = self.n_hba
            n_hba_2 = self.n_hba

        #create system
        subprocess.run('module load gromacs && gmx editconf -f {}.pdb -o {}_newbox.pdb -c -d {} -bt cubic'.format(local_pdb_file, local_pdb_file, buffer_size), shell=True)

        if n_hba_1 != 0 and n_hba_2 != 0:
            subprocess.run('module load gromacs && gmx insert-molecules -seed {} -f {}_newbox.pdb -ci {} -nmol {} -try 1000 -o out.pdb'.format(seed, local_pdb_file, self.hba_file_1, n_hba_1), shell=True)
            subprocess.run('module load gromacs && gmx insert-molecules -seed {} -f out.pdb -ci {} -nmol {} -try 1000 -o hba_out.pdb'.format(seed, self.hba_file_2, n_hba_2), shell=True)
        if n_hba_1 != 0 and n_hba_2 == 0:
            subprocess.run('module load gromacs && gmx insert-molecules -seed {} -f {}_newbox.pdb -ci {} -nmol {} -try 1000 -o hba_out.pdb'.format(seed, local_pdb_file, self.hba_file_1, n_hba_1), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        if n_hba_1 == 0 and n_hba_2 != 0:
            subprocess.run('module load gromacs && gmx insert-molecules -seed {} -f {}_newbox.pdb -ci {} -nmol {} -try 1000 -o hba_out.pdb'.format(seed, local_pdb_file, self.hba_file_2, n_hba_2), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        if self.n_hbd != 0:
            if n_hba_1 == 0 and n_hba_2 == 0:
                subprocess.run('module load gromacs && gmx insert-molecules -seed {} -f {}_newbox.pdb -ci {} -nmol {} -try 1000 -o hbd_out.pdb'.format(seed, local_pdb_file, self.hbd_file, self.n_hbd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
            else:
                subprocess.run('module load gromacs && gmx insert-molecules -seed {} -f hba_out.pdb -ci {} -nmol {} -try 1000 -o hbd_out.pdb'.format(seed, self.hbd_file, self.n_hbd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        if self.n_h2o != 0:
            if self.n_hbd == 0:
                subprocess.run('module load gromacs && gmx solvate -cp hba_out.pdb -maxsol {} -o protein_des.pdb'.format(self.n_h2o), shell=True)
            else:
                subprocess.run('module load gromacs && gmx solvate -cp hbd_out.pdb -maxsol {} -o protein_des.pdb'.format(self.n_h2o), shell=True)
        if self.n_h2o == 0:
            if self.n_hbd == 0:
                subprocess.run('mv hba_out.pdb protein_des.pdb'.format(self.n_h2o), shell=True)
            else:
                subprocess.run('mv hbd_out.pdb protein_des.pdb'.format(self.n_h2o), shell=True)
            
        
        # subprocess.run('rm out.pdb out2.pdb out3.pdb', shell=True)
        
    def apply_forcefield(self):
        """
        Apply the forcefield to the system

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        def remove_position_restraints(filename):
            with open(filename, "r") as f:
                lines = f.readlines()
                f.close()
            with open(filename, "w") as f:
                for line in lines[:-5]:
                        f.write(line)
                f.close()

        os.chdir(self.data_dir)
        shutil.copytree('../../data/forcefield/charmm36-jul2022.ff', './charmm36-jul2022.ff')

        subprocess.run('module load gromacs && gmx pdb2gmx -f protein_des.pdb -ff charmm36-jul2022 -water tip3p -o protein_des.gro', shell=True)
        if os.path.exists('./topol_Other2.itp'):
            remove_position_restraints("topol_Other2.itp")
        if os.path.exists('./topol_Other_chain_C.itp'):
            remove_position_restraints("topol_Other_chain_C.itp")
        if os.path.exists('./topol_Other.itp'):
            remove_position_restraints("topol_Other.itp")

    def run_sim(self, em_mdp, nvt_mdp, npt_mdp, production_mdp):
        """
        Placeholder function to show example docstring (NumPy format).

        Replace this function and doc string for your own project.

        Parameters
        ----------
        with_attribution : bool, Optional, default: True
            Set whether or not to display who the quote is from.

        Returns
        -------
        quote : str
            Compiled string including quote and optional attribution.
        """
        def run_em():
            subprocess.run('module load gromacs && gmx grompp -f {} -c protein_des.gro -p topol.top -o em.tpr'.format(em_mdp), shell=True)
            subprocess.run('module load gromacs && gmx mdrun -v -deffnm em -nt 0', shell=True)
        def run_nvt():
            subprocess.run('module load gromacs && gmx grompp -f {} -c em.gro -r em.gro -p topol.top -o nvt.tpr'.format(nvt_mdp), shell=True)
            subprocess.run('module load gromacs && gmx mdrun -v -deffnm nvt -nt 0', shell=True)
        def run_npt():
            subprocess.run('module load gromacs && gmx grompp -f {} -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr'.format(npt_mdp), shell=True)
            subprocess.run('module load gromacs && gmx mdrun -v -deffnm npt -nt 0', shell=True)
        def run_production():
            subprocess.run('module load gromacs && gmx grompp -f {} -c npt.gro -p topol.top -o production.tpr'.format(production_mdp), shell=True)
            subprocess.run('module load gromacs && gmx mdrun -v -deffnm production -nt 0', shell=True)

        run_em()
        run_nvt()
        run_npt()
        run_production()

    def analyze_trajectory(self, descriptor = 'gyrate'):
        """
        A function to analyize the trajectory and figure out the radius of gyration and possibliy other descriptors

        Parameters
        ----------
        with_attribution : bool, Optional, default: True
            Set whether or not to display who the quote is from.

        Returns
        -------
        RoG: float
            Radius of gyration of the protein
        """
        if not os.path.exists('./production.xtc'):
            raise Exception('The production run was unsuccessful, the .xtc file does not exist')

        # subprocess.run('module load gromacs && gmx energy -s production.tpr -f production.xtc -o gyrate.xvg', shell=True)
        if descriptor == 'gyrate':
            subprocess.run('module load gromacs && echo echo "1 0 " | gmx gyrate -s production.tpr -f production.xtc -o gyrate.xvg', shell=True)
            fname = 'gyrate.xvg'
        if descriptor == 'rmsd':
            subprocess.run('module load gromacs && echo echo "1 0 " | gmx rms -s npt.gro -f production.xtc -o rmsd.xvg', shell=True)
            fname = 'rmsd.xvg'
        with open(fname) as f:
            lines = f.readlines()
            f.close()
        
        values = []
        start = 0
        for i, l in enumerate(lines):
            if "s3 legend" in l:
                start = i
                break
        for l in lines[start+1:]:
            values.append(list(map(float, l.split()[1:]))[1])
        
        return np.mean(np.array(values))
        # running = df['Rg'].rolling(100).mean()

if __name__ == "__main__":
    pass
