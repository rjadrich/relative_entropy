max_loop=1000000 #set large as script terminates from within loop via errors or the subordinate max steps in the RE code
counter=1 #counter for loop

chmod 744 ./multi_g_rdf #MAKE SURE THIS IS EXECUTABLE RIGHT OFF THE BAT

while [ $counter -lt $max_loop ]
do


   ##############################################
   #########RELATIVE ENTROPY CALCULATION#########
   ##############################################

   ./relative_entropy 
   re_exit=$?

   if [[ $re_exit != 0 ]]
   then 
      echo "killing program"
      exit $re_exit
   fi

   ##############################################
   ######CALL SIM PROGRAM AND CALCULATE RDF######
   ##############################################

   chmod 744 ./auxillary.sh
   ./auxillary.sh
   aux_exit=$?

   if [[ $aux_exit != 0 ]]
   then
      echo "killing program"
      exit $aux_exit
   fi


   #read -p "Press any key to continue... " -n1 -s

done

