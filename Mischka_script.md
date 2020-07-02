inbreeding coefficient and multidimensional scaling (MDS) calculation
```bash
plink -bfile $input.bed –ibc -out ibc.out
plink -bfile $input.bed –mds-plot -out ibc.out
```

generate the gvcf file for each CH and ISR individual
```bash
gatk HaplotypeCaller \
	-I $bam_file \
	-O $sample_id.$chr.gvcf \
	-R $ref.fa \
	-L $chr \
	-ERC GVCF
```

combined the indivudual GVCF files 
```bash
gatk CombineGVCFs \
-R $genome.fa \
-O CombineGVCFs/all_dogs.g.vcf.gz
-V $dog1.gvcf \
-V $dog2.gvcf \
....
```

genotype the combined GVCFs
```bash
gatk GenotypeGVCFs \
-R /crex/proj/uppstore2017228/KLT.02.CCAN/b2017226_nobackup/hic_mischka/cf4.b6.12_analysis/chromium_ref/refdata-cf4.b6.12.IDs2.10XRef/fasta/genome.fa \
-V all_dogs.g.vcf.gz \
-O all_dogs.vcf.gz
```

Variants QC filtering
```bash
gatk SelectVariants \
-R /crex/proj/uppstore2017228/KLT.03.CGEN/Mischka/cf4_assemblies/cf4.b6.12.shortID.fasta \
-V $input_file \
-O $output_name.SNP.filtered.vcf \
-select "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--invertSelect \
-select-type SNP \
--restrict-alleles-to BIALLELIC


gatk SelectVariants \
-R /crex/proj/uppstore2017228/KLT.03.CGEN/Mischka/cf4_assemblies/cf4.b6.12.shortID.fasta \
-V $input_file \
-O $output_name.indel.filtered.vcf \
-select "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--invertSelect \
-select-type INDEL \
--restrict-alleles-to BIALLELIC
```

falcon denovo assemlby configs
```bash
[General]
input_fofn = input.fofn

length_cutoff = 8000
length_cutoff_pr = 14000

target = assembly

pa_DBsplit_option = -x500 -s400
ovlp_DBsplit_option = -x500 -s400

falcon_sense_option = --output_multi --min_idt 0.70 --min_cov 4 --max_n_read 200 --n_core 8 --min_cov_aln 4 --min_len_aln 40

overlap_filtering_setting = --max_diff 110 --max_cov 165 --min_cov 3 --n_core 8 --bestn 10

pa_HPCdaligner_option =   -v -B128 -t16 -e.70 -l1000 -s1000 -M28
ovlp_HPCdaligner_option = -v -B128 -t32 -e.96 -l500  -s1000 -M28 -h60

falcon_sense_skip_contained = false
skip_checks = True
dust = false
dazcon = false
```

ARCS+Links scaffolding with 10X linked reads
```bash
~/arcs/Examples/arcs-make arcs draft=$assembly reads=$barcoded.fq.gz m=20-10000 t=16
LINKS -f $assembly.fa -s empty.fof -a 0.9 -b assembly_c5_m20-10000_s98_r0.05_e30000_z500_l5_a0.9
```

repetitive elements annotation with EpeatMasker
```bash
RepeatMasker -pa 16 -s -species dog -dir ./out $genome.fasta
```

compare the assembly sequence with a linkage map 
```bash
chromonomer -p $linkage_map.tsv \
    -o ~/out -s $markers.sam \
    -a $final.assembly.agp
```

ARROWb long reads polishing

```bash
./variantCaller --algorithm=arrow -j 17  $pacbio_subreads.sort.bam.xml \
	-r $GSD1.0.fasta.xml \
	-o $variants.gff \
	-o $consensus.fastq \
	-o $consensus.fasta
```

blat map the marker sequence to assembly
```bash
blat $GSD1.0.fasta $markers.fa $position.psl
```

Juicer pipeline for re-map the HiC reads to the HiRise assembly
```bash
./juicer.sh \
-z $Hirise_out.fa \
-y $Hirise_out.fa_MboI.txt \
-D $juicer_run_path/ \
-d $juicer_run_path/ \
-p $Hirise_out.fa.genomeSize \
-q core -l core
```

Busco evaluate the completeness of genome assembly
```bash
run_BUSCO.py -c 5 -i $GSD1.0.fa -o $result.output -l $BUSCO_LINEAGE_SETS/mammalia_odb9 -m geno
```

BWA reference index and Mem mapping
```bash
##index the genome fasta file
bwa index $GSD1.0.fa
##mapping the pairend-ends reads in reference
bwa mem -t10 $GSD1.0.fa $reads_001_R1.fastq.gz $reads_001_R2.fastq.gz > $outputN.sam
```

PBjelly for gap filling
```xml
<jellyProtocol>
    <reference>$assembly.fasta</reference>
    <outputDir>$output_dir/</outputDir>
    <blasr>--minMatch 8 --sdpTupleSize 8 --minPctIdentity 75 --bestn 1 --nCandidates 10 --maxScore -500 --nproc 16 --noSplitSubreads</blasr>
    <input baseDir="$pb_subreads_fastq/">
        <job>$pb_subreads.fastq</job>
    </input>
</jellyProtocol>
```

```bash
Jelly.py setup protocol.xml
Jelly.py support protocol.xml
Jelly.py extraction protocol.xml
Jelly.py assembly protocol.xml -x "--nproc=16"
Jelly.py output protocol.xml
```

Pilon polishing with short reads
```bash
java -Xmx100G -jar $PILON_HOME/pilon.jar \
	--genome $assembly.fa \
	--frags $reads.10x.bam \
	--output $assembly_polished \
	--outdir $out_polishing/ \
	--fix bases \
	--changes \
	--diploid \
	--threads 20
```

Mappability with different Kmers GEM
```bash
##index the genome fasta
gem-indexer -T 16 -i $GSD1.0.fa -o $GSD1.0.gem_index
##calculate the mappability,$size is the Kmer size
gem-mappability -T 16 -I $GSD1.0.gem_index.gem -l $size -o $GSD1.0.gem_index.gem.$size
```

bamCoverage
```bash
bamCoverage -b $bam_file -o $sample_id.dep.bedGraph -p 4 --outFileFormat bedgraph --samFlagExclude 256 --binSize 25

bamCoverage -b $bam_file -o $sample_id.dep.bedGraph -p 4 --outFileFormat bedgraph --samFlagExclude 256 --binSize 25 --minMappingQuality 10
```