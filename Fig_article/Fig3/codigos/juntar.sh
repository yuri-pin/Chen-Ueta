for i in {1..40}; do
   let "j=$i-1" 
   cat  bif_x_-$j.dat bif_x_$i.dat > bif_x_-$i.dat
   echo 
   rm bif_x_$i.dat
   rm bif_x_-$j.dat
done
