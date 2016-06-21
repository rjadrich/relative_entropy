grompp -f grompp.mdp -po md_out.mdp -c conf.gro -n index.ndx -p topol.top -o topol.tpr

grompp_exit=$?

exit $grompp_exit
