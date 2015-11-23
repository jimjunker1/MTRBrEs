#to run type bash sh/cleanCSV.sh
cd data/
for file in *.csv
do
   perl -pi -e 's/\r/\n/g' $file
done
cd ..
