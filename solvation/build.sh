gmx editconf -f 1OIL_A_no_CAL.pdb -o 1OIL_newbox.pdb -c -d 2.0 -bt cubic
gmx insert-molecules -f 1OIL_newbox.pdb -ci chol.pdb -nmol 504 -o out.pdb
gmx insert-molecules -f out.pdb -ci cla.pdb -nmol 500 -o out2.pdb
gmx insert-molecules -f out2.pdb -ci urea.pdb -nmol 1000 -o protein_des.pdb
rm out.pdb out2.pdb 1OIL_newbox.pdb

gmx pdb2gmx -f protein_des.pdb -o protein_des.gro

gmx grompp -f em.mdp -c protein_des.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

#gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
#gmx mdrun -v -deffnm nvt

#gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
#gmx mdrun -deffnm npt

#gmx grompp -f production.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
#gmx mdrun -v -deffnm md_0_1
