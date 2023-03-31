Query = '"Betaherpesvirinae"[Organism] OR Betaherpesvirinae[All Fields]'
Command = f'esearch -db Nucleotide -query {Query} | efetch -format fasta > BetaherpesvirinaeData.txt'
#The command will use esearch to read the query. The -db will specfiy the datatbase. Query will connect the search with esearch.
#efetch will also work with esearch. The format will then format the file into an output file, which will be BetaherpesvirinaeData.txt

os.system(Command)

InFile = 'BetData.txt'
OutFile = 'Betaherpesvirinae'
database_type = 'nucl'
title = 'Betaherpesvirinae'

os.system(f'makeblastdb -in {input_file} -out {output_name} -title {title} -dbtype {database_type}')
# In will be for the input file. Out will be for the output file. Title will be the name of the database while dbtyep of fpr the type that was used. 

#Blasting
#Uses longest contig
inputFile = 'Longest.txt'
outputFile = 'BetResults.csv'
Blast = 'blastn -query ' + inputFile + ' -db Betaherpesvirinae -out ' + outputFile + ' -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitles" > blastout.tsv '
print(Blast)

os.system(Blast)

#making table of blast command outputs
os.system('echo "sacc pident length qstart qend sstart send bitscore evalue stitle" | cat - BetResults.csv > BetResults.tsv')
os.system('head -n 11 BetResults.tsv  > BetResults2.tsv')
