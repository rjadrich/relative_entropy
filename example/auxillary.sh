cd ./step_500

echo "//////////////////BEGIN AUX SCRIPT/////////////////" >> ../log.txt
echo "preparing simulation for step 500" >> ../log.txt

##############################################
################GROMMP SCRIPT#################
##############################################

echo "executing grompp to prepare for gromacs run" >> ../log.txt
chmod 744 ./grompp_files.sh
./grompp_files.sh

grompp_exit=$?

if [[ $grompp_exit != 0 ]]
then
   echo "grompp was unsuccessful -> killing!" >> ../log.txt
   echo "///////////////////END AUX SCRIPT//////////////////" >> ../log.txt
   echo " " >> ../log.txt
   exit $grompp_exit
fi

##############################################
################MDRUN SCRIPT##################
##############################################

echo "executing gromacs run" >> ../log.txt
chmod 744 ./mdrun_gromacs.sh
./mdrun_gromacs.sh

mdrun_exit=$?

if [[ $mdrun_exit != 0 ]]
then
   echo "mdrun was unsuccessful -> killing!" >> ../log.txt
   echo "///////////////////END AUX SCRIPT//////////////////" >> ../log.txt
   echo " " >> ../log.txt
   exit $mdrun_exit
fi

##############################################
#################RDF SCRIPT###################
##############################################

echo "calculating rdf with gromacs" >> ../log.txt
chmod 744 ./rdf_gromacs.sh
./rdf_gromacs.sh

rdf_exit=$?

if [[ $rdf_exit != 0 ]]
then
   echo "rdf calculation was unsuccessful -> killing!" >> ../log.txt
   echo "///////////////////END AUX SCRIPT//////////////////" >> ../log.txt
   echo " " >> ../log.txt
   exit $rdf_exit
fi

##############################################
#################WRITE DONE###################
##############################################

echo "creating done file" >> ../log.txt
echo " " >> ./done.txt

wd_exit=$?

if [[ $wd_exit != 0 ]]
then
   echo "done file write unsuccessful -> killing!" >> ../log.txt
   echo "///////////////////END AUX SCRIPT//////////////////" >> ../log.txt
   echo " " >> ../log.txt
   exit $wd_exit
fi

echo "///////////////////END AUX SCRIPT//////////////////" >> ../log.txt
echo " " >> ../log.txt

cd ..

