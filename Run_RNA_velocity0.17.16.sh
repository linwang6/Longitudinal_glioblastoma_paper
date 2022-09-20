# velocyto, version 0.17.16
# 

velocyto run -e SF112 --samtools-memory 150000 --samtools-threads 30 -b barcodes.tsv -m /refs/Homo_sapiens/velocity/hg38_repeats_repeatMasker_a
llTracks.gtf -o SF112_velocyto SF112.bam /refs/Homo_sapiens/velocity/genes.gtf
