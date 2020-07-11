
## get genomes

```
wget https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/\?HistoryId\=NCID_1_163497961_130.14.18.97_5555_1588507183_4102358427_0MetA0_S_HStore\&QueryKey\=1\&ReleaseType\=RefSeq\&FileType\=GENOME_FASTA\&Flat\=true
```

## make mmseqs indexes

```
mmseqs createdb ncbi-genomes-2020-05-03/*fna.gz seqDB
mmseqs createdb reverse_transcriptases/*fasta queries/reverse_transcriptases\nmmseqs createdb crispr_associated_proteins/*fasta queries/crispr_associated_proteins\nmmseqs createdb recombinases/*fasta queries/recombinases\nmmseqs createdb holliday_junction_resolvases/*fasta queries/holliday_junction_resolvases\nmmseqs createdb r2dm/r2dm.orf.fasta queries/r2dm\n
mmseqs createdb ncbi-genomes-2020-05-03/*fna.gz seqDB
mmseqs createdb crispr_associated_proteins/*fasta reverse_transcriptases/*fasta holliday_junction_resolvases/*fasta recombinases/*fasta r2dm/r2dm.orf.fasta queries/all
mmseqs createdb seqs/*fasta crisprDB
```

## run mmseqs search

```
mmseqs search --threads 48 -s 10 -a 1 --max-seqs 1000000 queries/all prokaryotes/seqDB alignments/all tmp
mmseqs search --threads 48 -s 10 -a 1 --max-seqs 1000000 CRISPRclass19/crisprDB prokaryotes/seqDB alignments/crispr tmp
```

These take a long time. Second one took a day.

## extract tables

Bed-like tables.

```
mmseqs convertalis --format-output theader,tstart,tend,qheader,evalue,pident queries/all prokaryotes/seqDB alignments/all alignments/all.bed.like
mmseqs convertalis --format-output target,tstart,tend,qheader,evalue,pident CRISPRclass19/crisprDB prokaryotes/seqDB alignments/crispr alignments/crispr.bed.like
```

Full tables for later inspection.

```
mmseqs convertalis --format-output query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,cigar,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov,qset,qsetid,tset,tsetid queries/all prokaryotes/seqDB alignments/all alignments/all.tsv
mmseqs convertalis --format-output query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,cigar,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov,qset,qsetid,tset,tsetid CRISPRclass19/crisprDB prokaryotes/seqDB alignments/crispr alignments/crispr.tsv
```

## extract BED files

CRISPR stuff

```
for f in cas{1..14}; do grep -i $f alignments/crispr.bed.like | tr ' ' '_' | cut -f 1-5 |awk '{ if ($2 > $3) { print $1, $3, $2, $4, $5 } else { print $0 }}' |tr ' ' '\t' | cut -f 1-5 | bedtools sort -i /dev/stdin >bed/$f.bed &; done
for f in cas1; do grep -i $f alignments/crispr.bed.like | grep -i -v cas10 | grep -i -v cas11 | grep -i -v cas12 | grep -i -v cas13 | grep -i -v cas14 | tr ' ' '_' | cut -f 1-5 |awk '{ i
f ($2 > $3) { print $1, $3, $2, $4, $5 } else { print $0 }}' |tr ' ' '\t' | cut -f 1-5 | bedtools sort -i /dev/stdin >bed/$f.bed &; done
```

Other stuff

```
grep -Ff <(cat reverse_transcriptases/*fasta | grep '^>' | cut -f 1 -d\| | tr -d '>') alignments/all.bed.like | tr ' ' '_'| cut -f 1-5 |awk '{ if ($2 > $3) { print $1, $3, $2, $4, $5 } else { print $0 }}' |tr ' ' '\t' | cut -f 1-5| bedtools sort -i /dev/stdin >bed/reverse_transcriptases.bed
grep r2dm alignments/all.bed.like | tr ' ' '_' | cut -f 1-5 |awk '{ if ($2 > $3) { print $1, $3, $2, $4, $5 } else { print $0 }}' |tr ' ' '\t' | cut -f 1-5| bedtools sort -i /dev/stdin >bed/r2dm.bed
grep -Ff <(cat holliday_junction_resolvases/*fasta | grep '^>' | cut -f 1 -d\| | tr -d '>') alignments/all.bed.like | tr ' ' '_'| cut -f 1-5 |awk '{ if ($2 > $3) { print $1, $3, $2, $4, 
$5 } else { print $0 }}' |tr ' ' '\t' | cut -f 1-5| bedtools sort -i /dev/stdin >bed/holliday_junction_resolvases.bed
grep -Ff <(cat recombinases/*fasta | grep '^>' | cut -f 1 -d\| | tr -d '>') alignments/all.bed.like | tr ' ' '_'| cut -f 1-5 |awk '{ if ($2 > $3) { print $1, $3, $2, $4, $5 } else { prin
t $0 }}' |tr ' ' '\t' | cut -f 1-5| bedtools sort -i /dev/stdin >bed/recombinases.bed
```

## find neighboring hits from different classes, filtering lightly for e-value

```
bedtools closest -d -a bed/reverse_transcriptases.bed -b bed/cas14.bed  | awk '$11 > 0 && $11 < 1000 && $5 < 1e-5 && $10 < 1e-5 { print }' | awk '{ if ($1 != f || ($3 < s || $2 > e)) { print; } f=$1; s=$2; e=$3; }' >closest/rt_cas14.bed
bedtools closest -d -a bed/r2dm.bed -b bed/cas14.bed  | awk '$11 > 0 && $11 < 1000 && $5 < 1e-5 && $10 < 1e-5 { print }' | awk '{ if ($1 != f || ($3 < s || $2 > e)) { print; } f=$1; s=$2
; e=$3; }' >closest/r2dm_cas14.bed
```

## update for cyanobacteria

The same steps above were repeated for cyanobacteria genomes from refseq, due to the fact that many cas14/RT hits were in that taxa.
This resulted in a second set of matches for each protein class.
The full set was then used to extract elements.

## element extraction

Here, we take candidate hits and extract ~10kb regions around them, then pull these into a pangenome graph and use that graph structure to help determine the boundaries of each element class.

Get the regions:

```
bedtools closest -d -a bed/reverse_transcriptases.bed -b bed/cas14.bed | awk '$11 > 0 && $11 < 2000 && $5 < 1e-10 && $10 < 1e-10 { print }' | awk '{ if ($1 != f || ($3 < s || $2 > e)) { print; } f=$1; s=$2; e=$3; }' | grep '^NC_\|^NZ_' | awk '{ start=$2; if ($7 < start) { start=$7 }; end=$3; if ($8 > end) { end=$8 }; l=end-start; if (end-start < 10000) { inc=int((10000-(end-start))/2); start-=inc; end+=inc;} if (start < 0) { end-=start; start = 0 }; print $1, start, end, $4, $5, $9, $10 }'  | tr ' ' '\t'  >regions/representative_rt_cas14_10kb.1.bed
bedtools closest -d -a bed/cyanobacteria_reverse_transcriptases.bed -b bed/cyanobacteria_cas14.bed | awk '$11 > 0 && $11 < 2000 && $5 < 1e-10 && $10 < 1e-10 { print }' | awk '{ if ($1 != f || ($3 < s || $2 > e)) { print; } f=$1; s=$2; e=$3; }' | grep '^NC_\|^NZ_' | awk '{ start=$2; if ($7 < start) { start=$7 }; end=$3; if ($8 > end) { end=$8 }; l=end-start; if (end-start < 10000) { inc=int((10000-(end-start))/2); start-=inc; end+=inc;} if (start < 0) { end-=start; start = 0 }; print $1, start, end, $4, $5, $9, $10 }'  | tr ' ' '\t'  >regions/cyanobacteria_rt_cas14_10kb.1.bed 
(cat regions/representative_rt_cas14_10kb.1.bed| grep -v -Ff <(cut -f 1 regions/cyanobacteria_rt_cas14_10kb.1.bed)  ; cat regions/cyanobacteria_rt_cas14_10kb.1.bed ) >regions/rt_cas14_10kb.bed
cat regions/rt_cas14_10kb.bed| parallel 'seq=$(echo {} | cut -f 1-2 -d _); target=$(echo {} | cut -f 2- ); echo $seq  $target' | tr ' ' '\t'  >regions/rt_cas14_10kb.clean.bed
```

Collect the sequences that we've hit.

First we make a sequence list:

```
cd prokaryotes
ls ncbi-genomes-2020-05-03/*.fna.gz cyanobacteria/GCF_*/*.fna.gz | parallel 'echo {} $(zcat {} | grep "^>" | tr -d ">" )' >seq_list.txt
```

Now, we grep out the files that contain hits and concatenate them into a master FASTA file:

```
cat regions/rt_cas14_10kb.bed| parallel 'seq=$(echo {} | cut -f 1-2 -d _); echo prokaryotes/$(grep $seq prokaryotes/seq_list.txt | head -1 | cut -f 1 -d \ )' | sort | uniq | parallel '
zcat {}' >genomes/merged/rt_cas14_10kb.merged.fa
samtools faidx genomes/merged/rt_cas14_10kb.merged.fa
```

And we query it for the regions we've selected, attaching the genus name to the sequence name as a prefix to aid in downstream interpretation of the data.

```
samtools faidx -n 1000000 genomes/merged/rt_cas14_10kb.merged.fa $(cat regions/rt_cas14_10kb.clean.bed | awk '{ print $1":"$2"-"$3 }' ) >regions/rt_cas14_10kb.clean.fa
paste <(cat regions/rt_cas14_10kb.bed ) <(cat regions/rt_cas14_10kb.clean.fa |tr -d '>' | paste - - | less -S ) | parallel -k 'genus=$(echo {} | cut -f 1 | cut -f 3 -d _); echo $genus {}' | awk '{ print ">"$1"."$9"_rt="$6"_cas14="$8; print $10 }' >regions/rt_cas14_10kb.clean.genus.fa
```

These regions can serve as the basis for a pangenome graph.
We drop short alignments <3kb to make sure that the homologies represented here are likely due to long matches indicative of mobile element families.

```
minimap2 -t 48 -c -w 1 -k 11 -X rt_cas14_10kb.clean.genus.fa rt_cas14_10kb.clean.genus.fa >rt_cas14_10kb.clean.genus.paf
<rt_cas14_10kb.clean.genus.paf fpa drop -l 3000 >rt_cas14_10kb.clean.genus.drop3k.paf
seqwish -s rt_cas14_10kb.clean.genus.fa -p rt_cas14_10kb.clean.genus.drop3k.paf -g rt_cas14_10kb.clean.genus.drop3k.gfa -t 48
```

Now, we manipulate the graph to remove tips that aren't supported by at least three sequences.
This will remove the flanking sequences from all of our repeats, resulting in paths embedded in the graph that match the minimal sequence of each element family.

```
odgi build -g rt_cas14_10kb.clean.genus.drop3k.gfa -o - \
    | odgi sort -i - -o - -p bSn -t 48 -O 10000 -Q 10000 -A >rt_cas14_10kb.clean.genus.drop3k.odgi
odgi viz -i rt_cas14_10kb.clean.genus.drop3k.odgi -o rt_cas14_10kb.clean.genus.drop3k.odgi.png -x 5000 -y 500 -R
```

Then, prune and remove exact duplicates which resulted from close overlapping regions in the input finding the same element sequence.
(Exact duplicates are possible, but unlikely as a few changes might be expected between insertions. Here the concern is to find a collection of sequences so we remove the dups.)

```
odgi prune -i rt_cas14_10kb.clean.genus.drop3k.odgi -o - -T -m 3 \
    | odgi paths -f -i - | paste - - \
    | awk '{ print length($2), $1, $2 }' \
    | sort -n | uniq -f 2 \
    | awk '$1 > 2000 { print $2; print $3 }' >rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.fa
```

This build/prune process can be iterated several more times to get a cleaner set of sequences in the graph:

```
minimap2 -c -w 1 -k 11 -X rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.fa rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.fa >rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.paf
seqwish -s rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.fa -p <(fpa drop -l 3000 <rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.paf ) -g rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.gfa
odgi build -g rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.gfa -o - \
    | odgi prune -i - -o - -T -m 3 \
    | odgi paths -i - -f  | paste - - \
    | awk '{ print length($2), $1, $2 }' \
    | sort -n | uniq -f 2 \
    | awk '$1 > 2000 { print $2; print $3 }' >rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.fa
minimap2 -c -w 1 -k 11 -X rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.fa rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.fa >rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.paf
seqwish -s rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.fa -p <(fpa drop -l 3000 <rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.paf ) -g rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.gfa
odgi build -g rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.gfa -o - | odgi sort -i - -o - -p ebn -A >rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi
odgi viz -i rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi -o rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi.png -x 4000 -y 500 -R
```

I've based my analysis on the final graph.

```
odgi paths -i rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi -H -D . | sed 's/path.length/path.length\trt.e.val\tcas14.e.val/' | tr -d '[' | tr -d ']' | awk 'NR == 1 { print } NR > 1 { $2 = $1"."$2; print }' | tr ' ' '\t' | sed 's/_rt=/\t/' | sed 's/_cas14=/\t/' | gzip >rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi.paths.tsv.gz
odgi bin -i rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi -w 50 -D . | sed 's/path.name/path.name\trt.e.val\tcas14.e.val/' | tr -d '[' | tr -d ']' | tr ' ' '\t' | sed 's/_rt=/\t/' | sed 's/_cas14=/\t/' | sed 's/\t[0-9\.]*E-/\t/g' | gzip >rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi.bin50.tsv.gz
```

In R, these can be examined:

```R
require(tidyverse)
require(ape)      
require(phyclust) 
require(ggfortify)
require(ggtree)

casRTs <- read.delim('rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi.paths.tsv.gz')
casRTs.matrix <- casRTs[ , !names(casRTs) %in% c("group.name","path.name","path.length","node.count","rt.e.val","cas14.e.val")]
casRTs.dist <- dist(casRTs.matrix)
casRTs.tree <- nj(casRTs.dist)

ggtree(casRTs.tree) %<+% data.frame(node=1:(nrow(casRTs.tree$edge)+1), path.name=factor(c(as.character(casRTs$path.name), rep("",nrow(casRTs.tree$edge)-nrow(casRTs)+1))), group.name=factor(c(as.character(casRTs$group.name),rep("internal",nrow(casRTs.tree$edge)-nrow(casRTs)+1)), levels=c(levels(casRTs$group.name),"internal"))) + aes(color=group.name, label=path.name) + scale_color_manual("genus", values=c(rainbow(20),"black")) + geom_text(nudge_x=5,size=3); ggsave("rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi.paths.tree.pdf", height=20, width=20)


# todo change this to the right file
casRT.bin50 <- read.delim('rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean+Cas14_genes.odgi.bin50.tsv.gz')
casRT.bin50$path.name <- fct_inorder(casRT.bin50$path.name)

ggplot(casRT.bin50, aes(x=bin, y=path.name, fill=path.prefix)) + geom_tile() + scale_fill_manual("genus", values=c(rainbow(19),"black")); ggsave('rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi.bin50.tile.pdf', height=20, width=20)
ggplot(casRT.bin50, aes(x=bin, y=path.name, color=cas14.e.val*rt.e.val, fill=cas14.e.val*rt.e.val)) + geom_tile(); ggsave("rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi.bin50.tile.cas14_x_rt_e-val.pdf", height=20, width=20)
ggplot(casRT.bin50, aes(x=bin, y=path.name, color=cas14.e.val+rt.e.val, fill=cas14.e.val+rt.e.val)) + geom_tile(); ggsave("rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi.bin50.tile.cas14_+_rt_e-val.pdf", height=20, width=20)
ggplot(casRT.bin50, aes(x=bin, y=path.name, color=cas14.e.val, fill=cas14.e.val)) + geom_tile(); ggsave("rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi.bin50.tile.cas14_e-val.pdf", height=20, width=20)
ggplot(casRT.bin50, aes(x=bin, y=path.name, color=rt.e.val, fill=rt.e.val)) + geom_tile(); ggsave("rt_cas14_10kb.clean.genus.drop3k.prune-T-m3.uniq.reclean.odgi.bin50.tile.rt_e-val.pdf", height=20, width=20)
```