qtlseq -o test \
       -n1 20 \
       -n2 20 \
       -w 100 \
       -s 20 \
       -r qtlseq_ref.fasta \
       -p qtlseq_parent.1.fastq.gz,qtlseq_parent.2.fastq.gz \
       -b1 qtlseq_bulk1.1.fastq.gz,qtlseq_bulk1.2.fastq.gz \
       -b2 qtlseq_bulk2.1.fastq.gz,qtlseq_bulk2.2.fastq.gz
