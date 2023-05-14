# CAN_cotton
Cold associated natural antisense transcripts (CAN)  were identified from G. hirsutum and G. barbadense in this project

Step 1. Building genome index

	extract_splice_sites.py TM-1_V2.1.gene.gtf > TM-1_V2.1.ss
	extract_exons.py TM-1_V2.1.gene.gtf > TM-1_V2.1.exon
	hisat2-build --ss TM-1_V2.1.ss --exon TM-1_V2.1.exon TM-1_V2.1.fa TM-1_V2.1
	
	extract_splice_sites.py Hai7124_V1.1.gene.gtf > Hai7124_V1.1.ss
  	extract_exons.py Hai7124_V1.1.gene.gtf > Hai7124_V1.1.exon
  	hisat2-build --ss Hai7124_V1.1.ss --exon Hai7124_V1.1.exon Hai7124_V1.1.fa Hai7124_V1.1

Step 2. Read aligning


	ls TM1*gz|perl -ne 's/.R[12].clean.fastq.gz//;print' | uniq |while read id; do hisat2 -p 8 --dta --rna-strandness RF -x ../../ref/TM-1_V2.1 -1 ${id}.R1.clean.fastq.gz -2 ${id}.R2.clean.fastq.gz -S ${id}.sam

	ls H7*gz|perl -ne 's/.R[12].clean.fastq.gz//;print' | uniq |while read id; do hisat2 -p 8 --dta --rna-strandness RF -x ../../ref/Hai7124_V1.1 -1 ${id}.R1.clean.fastq.gz -2 ${id}.R2.clean.fastq.gz -S ${id}.sam
	
	ls *sam | while read id; do samtools sort -@ 1 --output-fmt BAM -o ${id}.sorted.bam $id 

Step 3. assemble

	ls *bam | while read id; do stringtie --rf -f 0.1 -j 10 -c 10  -p 2 -o ${id}.gtf $id;done

Step 4.  lncRNA prediction

	######## annonation lncRNA #############
	
	ls TM1*sorted.bam.gtf > TM1_mergelist.txt
	ls H7*sorted.bam.gtf > H7_mergelist.txt
	
	stringtie --merge -o TM1_merged_1.gtf -c 10 -f 0.1  -G TM1.gene.gtf ./TM1_mergelist.txt
	stringtie --merge -o H7_merged_1.gtf -c 10 -f 0.1  -G H7.gene.gtf ./H7_mergelist.txt
	
	gffcompare -r ~/project/ref/TM-1_V2.1.gene.gtf  -p 4 TM1_merged_1.gtf -o TM1_merged_lncRNA
	gffcompare -r ~/project/ref/Hai7124_V1.1.gene.gtf  -p 4 H7_merged_1.gtf -o H7_merged_lncRNA
	
	
	awk '$3 == "i" || $3 == "y" || $3 == "u" {print $0}' TM1_merged_lncRNA.merged_1.gtf.tmap > TM1_novel.gtf.tmap
	awk '$10 >200 {print}' TM1_novel.gtf.tmap > TM1_novel.longRNA.gtf.tmap
	perl ~/bin/extract_gtf_by_name.pl TM1_novel.longRNA.gtf.tmap 4 TM1_merged_1.gtf  > TM1_novel.longRNA.gtf
	gffread -g ~/project/ref/TM-1_V2.1.fa -w TM1_exon.fa ./TM1_novel.longRNA.gtf
	TransDecoder.LongOrfs -t TM1_exon.fa # This step generated a file named longest_orfs.pep
	pfam_scan.pl -cpu 8 -fasta ./TM1_exon.fa.transdecoder_dir/longest_orfs.pep -dir ~/software/pfam/ > TM1_pfam_scan.txt
	perl -ne 'print if /Domain/' TM1_pfam_scan.txt |perl -ne '/^(\w+\.\w+)\.\w+\.\w+/;print "$1\n" unless exists $hash{$1};$hash{$1}++ ' > TM1_novel.transcript_with_domain.txt
	CPC2.py -i TM1_exon.fa -o TM1_cpc_output.txt
	perl -ne 'print if /noncoding/' TM1_cpc_output.txt |perl -ne '/^(\w+\.\w+)\.\w+/; print "$1\n" '  > TM1_novel.transcript_cpc_nocoding.txt
	cat TM1_novel.transcript_cpc_nocoding.txt TM1_novel.transcript_with_domain.txt |sort|uniq -d > TM1_intersection.txt
	sort TM1_novel.transcript_cpc_nocoding.txt TM1_intersection.txt  |uniq -u > TM1_lncRNA_list.txt
	perl ~/bin/extract_gtf_by_geneid.pl TM1_lncRNA_list.txt TM1_novel.longRNA.gtf  > TM1_lncRNA.gtf
	perl ~/bin/rename_lncRNA.pl TM1_lncRNA.gtf Gh > TM1_lncRNA.rename.gtf
	cat  TM1_lncRNA.rename.gtf ~/project/ref/TM-1_V2.1.gene.gtf > TM1_lncRNA_mRNA.gtf
	
	awk '$3 == "i" || $3 == "y" || $3 == "u" {print $0}' H7_merged_lncRNA.merged_1.gtf.tmap > H7_novel.gtf.tmap
	awk '$10 >200 {print}' H7_novel.gtf.tmap > H7_novel.longRNA.gtf.tmap
	perl ~/bin/extract_gtf_by_name.pl H7_novel.longRNA.gtf.tmap 4 H7_merged_1.gtf  > H7_novel.longRNA.gtf
	gffread -g ~/project/ref/Hai7124_V1.1.fa -w H7_exon.fa ./H7_novel.longRNA.gtf
	TransDecoder.LongOrfs -t H7_exon.fa # This step generated a file named longest_orfs.pep
	pfam_scan.pl -cpu 8 -fasta ./H7_exon.fa.transdecoder_dir/longest_orfs.pep -dir ~/software/pfam/ > H7_pfam_scan.txt
	perl -ne 'print if /Domain/' H7_pfam_scan.txt |perl -ne '/^(\w+\.\w+)\.\w+\.\w+/;print "$1\n" unless exists $hash{$1};$hash{$1}++ ' > H7_novel.transcript_with_domain.txt
	CPC2.py -i H7_exon.fa -o H7_cpc_output.txt
	perl -ne 'print if /noncoding/' H7_cpc_output.txt |perl -ne '/^(\w+\.\w+)\.\w+/; print "$1\n" '  > H7_novel.transcript_cpc_nocoding.txt
	cat H7_novel.transcript_cpc_nocoding.txt H7_novel.transcript_with_domain.txt |sort|uniq -d > H7_intersection.txt
	sort H7_novel.transcript_cpc_nocoding.txt H7_intersection.txt  |uniq -u > H7_lncRNA_list.txt
	perl ~/bin/extract_gtf_by_geneid.pl H7_lncRNA_list.txt H7_novel.longRNA.gtf  > H7_lncRNA.gtf
	perl ~/bin/rename_lncRNA.pl H7_lncRNA.gtf Gb > H7_lncRNA.rename.gtf
	cat  H7_lncRNA.rename.gtf ~/project/ref/Hai7124_V1.1.gene.gtf > H7_lncRNA_mRNA.gtf
		
	sortBed -i TM1_lncRNA_mRNA.gtf > 1 && mv 1 TM1_lncRNA_mRNA.gtf
	gffread -T TM1_lncRNA_mRNA.gtf > 1 && mv 1 TM1_lncRNA_mRNA.gtf
	perl -ne 'next unless /\ttranscript\t/;/gene_id "([\.\w]+)"/;print unless exists $hash{$1};$hash{$1}++' TM1_lncRNA_mRNA.gtf > TM1_lncRNA_mRNA.transcript.uniq.lst
	perl ~/bin/uniq.transcript.pl TM1_lncRNA_mRNA.transcript.uniq.lst TM1_lncRNA_mRNA.gtf > TM1_lncRNA_mRNA.uniq.gtf
	
	sortBed -i H7_lncRNA_mRNA.gtf > 1 && mv 1 H7_lncRNA_mRNA.gtf
	gffread -T H7_lncRNA_mRNA.gtf > 1 && mv 1 H7_lncRNA_mRNA.gtf
	perl -ne 'next unless /\ttranscript\t/;/gene_id "([\.\w]+)"/;print unless exists $hash{$1};$hash{$1}++' H7_lncRNA_mRNA.gtf > H7_lncRNA_mRNA.transcript.uniq.lst
	perl ~/bin/uniq.transcript.pl H7_lncRNA_mRNA.transcript.uniq.lst H7_lncRNA_mRNA.gtf > H7_lncRNA_mRNA.uniq.gtf
	
	rm *novel.transcript_cpc_nocoding.txt
	rm *novel.transcript_with_domain.txt
	rm *pfam_scan.txt
	rm *intersection.txt
	
	
Step 5. NAT prediction
	
	######## annonation lncRNA #############
		
	gffcompare -r TM1_lncRNA_mRNA.uniq.gtf -p 4 TM1_merged_1.gtf -o TM1_merged_NAT
	gffcompare -r H7_lncRNA_mRNA.uniq.gtf -p 4 H7_merged_1.gtf -o H7_merged_NAT
	
	awk '$3 == "x" {print $0}' TM1_merged_NAT.merged_1.gtf.tmap > TM1_novel.NAT.gtf.tmap
	awk '$3 == "x" {print $0}' H7_merged_NAT.merged_1.gtf.tmap > H7_novel.NAT.gtf.tmap
	
	perl ~/bin/extract_gtf_by_name.pl TM1_novel.NAT.gtf.tmap 4 TM1_merged_1.gtf > TM1_novel.NAT.gtf 
	perl ~/bin/extract_gtf_by_name.pl H7_novel.NAT.gtf.tmap 4 H7_merged_1.gtf > H7_novel.NAT.gtf 
			
	perl ~/bin/rename_lncRNA.pl TM1_novel.NAT.gtf Gh_NAT > TM1_lncRNA.NAT.rename.gtf
	perl ~/bin/rename_lncRNA.pl H7_novel.NAT.gtf Gb_NAT > H7_lncRNA.NAT.rename.gtf
	
	cat  TM1_lncRNA.NAT.rename.gtf TM1_lncRNA_mRNA.uniq.gtf > TM1_lncRNA_mRNA.uniq.NAT.gtf
	cat  H7_lncRNA.NAT.rename.gtf H7_lncRNA_mRNA.uniq.gtf > H7_lncRNA_mRNA.uniq.NAT.gtf
	
	sortBed -i TM1_lncRNA_mRNA.uniq.NAT.gtf > 1 && mv 1 TM1_lncRNA_mRNA.uniq.NAT.gtf
	gffread -T TM1_lncRNA_mRNA.uniq.NAT.gtf > 1 && mv 1 TM1_lncRNA_mRNA.uniq.NAT.gtf
	
	sortBed -i H7_lncRNA_mRNA.uniq.NAT.gtf > 1 && mv 1 H7_lncRNA_mRNA.uniq.NAT.gtf
	gffread -T H7_lncRNA_mRNA.uniq.NAT.gtf > 1 && mv 1 H7_lncRNA_mRNA.uniq.NAT.gtf
	
	perl -ne 'next unless /\ttranscript\t/;/gene_id "([\.\w]+)"/;print unless exists $hash{$1};$hash{$1}++' TM1_lncRNA_mRNA.uniq.NAT.gtf > TM1_lncRNA_mRNA.uniq.NAT.lst
	perl ~/bin/uniq.transcript.pl TM1_lncRNA_mRNA.uniq.NAT.lst TM1_lncRNA_mRNA.uniq.NAT.gtf >1 && mv 1 TM1_lncRNA_mRNA_NAT.gtf
	
	
	perl -ne 'next unless /\ttranscript\t/;/gene_id "([\.\w]+)"/;print unless exists $hash{$1};$hash{$1}++' H7_lncRNA_mRNA.uniq.NAT.gtf > H7_lncRNA_mRNA.uniq.NAT.lst
	perl ~/bin/uniq.transcript.pl H7_lncRNA_mRNA.uniq.NAT.lst H7_lncRNA_mRNA.uniq.NAT.gtf >1 && mv 1 H7_lncRNA_mRNA_NAT.gtf
	
Step 6. Building collinearity anchors

	mkdir synteny
	python -m jcvi.formats.gff bed --type=transcript --key=geneID TM1_lncRNA_mRNA_NAT.gtf -o synteny/Gh.bed
	python -m jcvi.formats.gff bed --type=transcript --key=geneID Hai7124_lncRNA_mRNA_NAT.gtf -o synteny/Gb.bed
	
	cd synteny
	perl -ne 's/(GB_\w+)_[keq]+/$1/;print' Gb.bed > 1 && mv 1 Gb.bed 
	perl -ne 's/(GH_\w+)_[keq]+/$1/;print' Gh.bed > 1 && mv 1 Gh.bed 
		
	gffread -g ~/project/ref/Hai7124_V1.1.fa -w Hai7124_lncRNA_mRNA_NAT.fa Hai7124_lncRNA_mRNA_NAT.gtf 
	gffread -g ~/project/ref/TM1_V2.1.fa -w TM1_lncRNA_mRNA_NAT.fa TM1_lncRNA_mRNA_NAT.gtf 
	
	perl -ne 's/(GH_\w+)_[eqk]+/$1/;s/(>\w+)\.\d+/$1/;print' TM1_lncRNA_mRNA_NAT.fa >1 && mv  1 Gh.cds
	perl -ne 's/(GB_\w+)_[eqk]+/$1/;s/(>\w+)\.\d+/$1/;print' Hai7124_lncRNA_mRNA_NAT.fa >1 && mv 1 Gb.cds
	
	perl -ne 'print unless /NAT/' Gh.bed > 1 && mv 1 Gh.bed
	perl -ne 'print unless /NAT/' Gb.bed > 1 && mv 1 Gb.bed
	
	perl remove_NAT_from_fa.pl Gh.cds > 1 && mv 1 Gh.cds
	perl remove_NAT_from_fa.pl Gb.cds > 1 && mv 1 Gb.cds
	
	python -m jcvi.compara.catalog ortholog Gh Gb --cscore 0.9 --cpus 5 --no_strip_names
	
Step 7. reciprocal best hits
	
	makeblastdb -in Gh.cds -dbtype nucl
	makeblastdb -in Gb.cds -dbtype nucl
	
	blastn -query Gb.cds -db Gh.cds -outfmt 6 -out Gb_Gh.blastn.m6 -num_threads 10 -max_target_seqs 3 -evalue 1e-5 
	blastn -query Gh.cds -db Gb.cds -outfmt 6 -out Gh_Gb.blastn.m6 -num_threads 10 -max_target_seqs 3 -evalue 1e-5
	
	perl extract_RB_hit.pl Gh_Gb.blastn.m6 Gb_Gh.blastn.m6 > Gh_Gb_RB_hit.lst
	
Step 8. Obtain Gh_Gb conserved NATs

	perl creat_UpSet_data.pl Gh.Gb.anchors.lst Gh_Gb_RB_hit.lst > synteny_RBblast_Upset.data.txt
	
	perl -lane 'print "$F[0]" if $F[1] == 1 && $F[2] == 1' synteny_RBblast_Upset.data.txt > Gh_Gb.total.homology.lst
	perl -ne '/(GH_\w+)_(GB_\w+)/;print "$1\t$2\n"' Gh_Gb.total.homology.lst > 1 && mv 1 Gh_Gb.total.homology.lst
	
	perl extract_homology_NAT.pl TM1_lncRNA_mRNA_NAT.gtf H7_lncRNA_mRNA_NAT.gtf Gh_Gb.total.homology.lst > Gh_Gb_homology_NAT.lst



