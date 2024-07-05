#grep -oP "http.*>" Index_Proteomes.txt | grep  "[0-9]/>$" | tr -d '>' > Protein_list_clean.txt
count=1
last=`wc -l Proteins_clean.txt `
for  f in `cat  Proteins_clean.txt`
do
rm index.html
wget -q --no-parent $f  # https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Archaea/UP000000242/
file=` grep -oP "\".*.fasta.gz\"" index.html | grep -v DNA | tr -d '"'` 
wget "$f"$file
nfile=` zcat $file | head -n 1 | awk -F'=' '{gsub(" OX","",$0); gsub(" ","_",$2); print $2"_"'$count' }'`
echo " working on $nfile , number $count of $last" 
mv $file $nfile.fasta.gz >>issues.log
 let count++
 sleep 15 # nothing so powerful as a random sleep, they wont think we are a  machine  
 done 
