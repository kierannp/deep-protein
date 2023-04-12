gmx insert-molecules -ci pdbs/chol.pdb -nmol 1000 -box 25 25 25 -o out.pdb
gmx insert-molecules -f out.pdb -ci pdbs/cla.pdb -nmol 1000 -o out2.pdb
gmx insert-molecules -f out2.pdb -ci pdbs/urea.pdb -nmol 2000 -o des.pdb
rm out.pdb out2.pdb

gmx pdb2gmx -f des.pdb
gmx grompp -f em.mdp -c conf.gro -p topol.top -o em.tpr
gmx mdrun -s em.tpr -v
gmx grompp -f npt.mdp -c conf.gro -o npt.tpr
