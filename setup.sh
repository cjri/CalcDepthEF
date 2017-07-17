ls Data/*sam > Input_files.in
while read line
  do 
  sort -nk1 $line > temp
  mv temp $line
done < Input_files.in
head -n1 Input_files.in > temp1
while read line
  do 
  head -n50 $line | cut -f1 | grep : | head -n3 | tail -n2 > col1
done < temp1
perl getline.pl

