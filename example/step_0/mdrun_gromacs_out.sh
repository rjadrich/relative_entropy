mdrun -s topol.tpr -c conf_out.gro -o traj.trr -x traj.xtc

mdrun_exit=$?

exit $mdrun_exit
