dragen  := "./result/bin/dragen-os"

test_segfault: test_data build_reference
    {{dragen}} -r dragmap --RGSM test --num-threads 2 -1 SRX882903_T2.fastq.gz > output.log

build_reference: test_data
    # mkdir dragmap
    {{dragen}} \
        --build-hash-table true \
        --ht-reference GRCh38_chr21.fa \
        --output-directory dragmap \
        --ht-num-threads 2

test_data:
    wget -nc https://raw.githubusercontent.com/nf-core/test-datasets/nascent/testdata/SRX882903_T2.fastq.gz
    wget -nc https://raw.githubusercontent.com/nf-core/test-datasets/nascent/reference/GRCh38_chr21.fa

issue_10:
    # https://nfcore.slack.com/archives/CGFUX04HZ/p1718268624579699?thread_ts=1718108323.285709&cid=CGFUX04HZ
    # https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R1_001.fastq.gz
    # https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R2_001.fastq.gz
    {{dragen}} \
        -r /home/tigem/h.poddar/short_reads_pipelines/ref_genome/GRCh37/hash/ \
        -1 /home/tigem/h.poddar/DRAGMAP/2A1_CGATGT_L001_R1_001.fastq.gz \
        -2 /home/tigem/h.poddar/DRAGMAP/2A1_CGATGT_L001_R2_001.fastq.gz \
        | samtools view -h -O BAM - > aligned.bam
