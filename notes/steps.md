
# get genomes

wget https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/\?HistoryId\=NCID_1_163497961_130.14.18.97_5555_1588507183_4102358427_0MetA0_S_HStore\&QueryKey\=1\&ReleaseType\=RefSeq\&FileType\=GENOME_FASTA\&Flat\=true

# make mmseqs indexes

mmseqs createdb ncbi-genomes-2020-05-03/*fna.gz seqDB
mmseqs createdb reverse_transcriptases/*fasta queries/reverse_transcriptases\nmmseqs createdb crispr_associated_proteins/*fasta queries/crispr_associated_proteins\nmmseqs createdb recombinases/*fasta queries/recombinases\nmmseqs createdb holliday_junction_resolvases/*fasta queries/holliday_junction_resolvases\nmmseqs createdb r2dm/r2dm.orf.fasta queries/r2dm\n
mmseqs createdb ncbi-genomes-2020-05-03/*fna.gz seqDB
mmseqs createdb crispr_associated_proteins/*fasta reverse_transcriptases/*fasta holliday_junction_resolvases/*fasta recombinases/*fasta r2dm/r2dm.orf.fasta queries/all
mmseqs createdb seqs/*fasta crisprDB

# run mmseqs search

mmseqs search --threads 48 -s 10 -a 1 --max-seqs 1000000 queries/all prokaryotes/seqDB alignments/all tmp
mmseqs search --threads 48 -s 10 -a 1 --max-seqs 1000000 CRISPRclass19/crisprDB prokaryotes/seqDB alignments/crispr tmp

# extract tables

mmseqs convertalis --format-output theader,tstart,tend,qheader,evalue,pident queries/all prokaryotes/seqDB alignments/all alignments/all.bed.like
mmseqs convertalis --format-output target,tstart,tend,qheader,evalue,pident CRISPRclass19/crisprDB prokaryotes/seqDB alignments/crispr alignments/crispr.bed.like

mmseqs convertalis --format-output query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,cigar,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov,qset,qsetid,tset,tsetid queries/all prokaryotes/seqDB alignments/all alignments/all.tsv
mmseqs convertalis --format-output query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,cigar,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov,qset,qsetid,tset,tsetid CRISPRclass19/crisprDB prokaryotes/seqDB alignments/crispr alignments/crispr.tsv

# extract BED files

for f in cas{1..14}; do grep -i $f alignments/crispr.bed.like | tr ' ' '_' | cut -f 1-5 |awk '{ if ($2 > $3) { print $1, $3, $2, $4, $5 } else { print $0 }}' |tr ' ' '\t' | cut -f 1-5 | bedtools sort -i /dev/stdin >bed/$f.bed &; done
for f in cas1; do grep -i $f alignments/crispr.bed.like | grep -i -v cas10 | grep -i -v cas11 | grep -i -v cas12 | grep -i -v cas13 | grep -i -v cas14 | tr ' ' '_' | cut -f 1-5 |awk '{ i
f ($2 > $3) { print $1, $3, $2, $4, $5 } else { print $0 }}' |tr ' ' '\t' | cut -f 1-5 | bedtools sort -i /dev/stdin >bed/$f.bed &; done

grep -Ff <(cat reverse_transcriptases/*fasta | grep '^>' | cut -f 1 -d\| | tr -d '>') alignments/all.bed.like | tr ' ' '_'| cut -f 1-5 |awk '{ if ($2 > $3) { print $1, $3, $2, $4, $5 } else { print $0 }}' |tr ' ' '\t' | cut -f 1-5| bedtools sort -i /dev/stdin >bed/reverse_transcriptases.bed
grep r2dm alignments/all.bed.like | tr ' ' '_' | cut -f 1-5 |awk '{ if ($2 > $3) { print $1, $3, $2, $4, $5 } else { print $0 }}' |tr ' ' '\t' | cut -f 1-5| bedtools sort -i /dev/stdin >bed/r2dm.bed
grep -Ff <(cat holliday_junction_resolvases/*fasta | grep '^>' | cut -f 1 -d\| | tr -d '>') alignments/all.bed.like | tr ' ' '_'| cut -f 1-5 |awk '{ if ($2 > $3) { print $1, $3, $2, $4, 
$5 } else { print $0 }}' |tr ' ' '\t' | cut -f 1-5| bedtools sort -i /dev/stdin >bed/holliday_junction_resolvases.bed
grep -Ff <(cat recombinases/*fasta | grep '^>' | cut -f 1 -d\| | tr -d '>') alignments/all.bed.like | tr ' ' '_'| cut -f 1-5 |awk '{ if ($2 > $3) { print $1, $3, $2, $4, $5 } else { prin
t $0 }}' |tr ' ' '\t' | cut -f 1-5| bedtools sort -i /dev/stdin >bed/recombinases.bed

# find neighboring hits from different classes, filtering for 

bedtools closest -d -a bed/reverse_transcriptases.bed -b bed/cas14.bed  | awk '$11 > 0 && $11 < 1000 && $5 < 1e-5 && $10 < 1e-5 { print }' | awk '{ if ($1 != f || ($3 < s || $2 > e)) { print; } f=$1; s=$2; e=$3; }' >closest/rt_cas14.bed
bedtools closest -d -a bed/r2dm.bed -b bed/cas14.bed  | awk '$11 > 0 && $11 < 1000 && $5 < 1e-5 && $10 < 1e-5 { print }' | awk '{ if ($1 != f || ($3 < s || $2 > e)) { print; } f=$1; s=$2
; e=$3; }' >closest/r2dm_cas14.bed