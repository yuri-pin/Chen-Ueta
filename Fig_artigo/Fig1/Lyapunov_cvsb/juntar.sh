for i in {1..500}; do
   let "j=$i-1" 
   cat  x_-$j.dat x_$i.dat > x_-$i.dat
   echo 
   rm x_$i.dat
   rm x_-$j.dat
done
