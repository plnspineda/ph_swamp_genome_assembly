Supplementary Notes 1

Processing of the Illumina paired-end short-reads data

Before using the short-reads for downstream analysis, we trimmed the Illumina adapters first using Trim Galore using the following parameters:

	trim_galore --quality 30 --length 110 --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 –	three_prime_clip_R2 10 --paired read_1.fastq read_2.fastq -o output_dir

The short reads were primarily used for estimation of genome size and heterozygosity score using GenomeScope. Additionally, the reads were also used for assessment of the assembly with Merqury. To generate k-mer count (21) of the Illumina short-reads, we run it with meryl using the following parameters:

	meryl k=21 count output swamp_buf.meryl $reads
	meryl union-sum output swamp_buf.meryl *reads.meryl

The meryl file was uploaded to http://qb.cshl.edu/genomescope/genomescope2.0/ and analysed with default parameters for genome size and heterozygosity estimation.

Contigs, scaffolds and final assembly were all assessed with Merqury.

	merqury.sh short-reads.meryl asm.fasta QV

Genome assembly

Contig assembly

We run our de novo contig genome assembly using HiFiasm with the command line:

	hifiasm -o output_dir -t 32 *fastq.gz

Removing duplications and junks

Firstly, we map the assembly with HiFi long reads with this argument using minimap2:

	minimap2 -xasm20 -I3g -t32 asm.fasta reads.fasta | gzip -c - > aln.paf.gz

Next, we followed purge_dups pipeline with the following commands:

	# calculate read depth histogram and base-level read depth
	pbcstat *paf.gz
	calcuts PB.stat > cutoffs 2>calcuts.log

	# split the assembly and do self-alignment
	split_fa asm.fasta > asm.fasta.split
	minimap2 -xasm5 -DP asm.fasta.split asm.fasta.split | gzip -c - > asm.split.self.paf.gz

	# identify regions and sequences to the haplobin
	purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > dups.bed 2> purge_dups.log

	# get purged primary and haplotig sequences from draft assembly
	get_seqs -e dups.bed asm.fasta

	filterbyname.sh in=asm.fasta out=purged_asm.fasta names=junk.list exclude

Scaffolding

We used the Arima Mapping Pipeline Guidelines to map the HiC short-reads data against the contigs and produced a BAM alignment file for downstream processing. Error correction option was not performed to avoid removal of segmental duplications resolved by HiFiasm. Then we run scaffolding as below:

	yahs –no-contig-ec asm_purged.fasta  asm_purged.bam

Visualising HiC contact maps

We used Juicebox to visualise HiC contact maps and to join scaffolds with strong signals. To generate a .hic and .assembly file as inputs to Juicebox, we run with the following commands:

	# index
	samtools faidx asm.fasta

	# generate .assembly file
	./yahs/juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp asm.fasta.fai > 	out_JBAT.log 2>&1

	# create file with scaffolds sizes
	cut -f 1,2 yahs.out_scaffolds_final.fa.fai > scaffolds_final.sizes

	# generate .hic file
	java -jar juicer_tools.1.8.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part < (cat 	out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}') && mv out_JBAT.hic.part 	out_JBAT.hic

The output file can then be loaded to Juicebox for visualisation and manual curation.

Alignment with reference genome assembly

We determine the chromosome number of the swamp buffalo by comparing homology with the riverine genome (UOA_WB_1). Alignment was done using meryl for counting k-mer of the assembly and mapping with winnowmap:

	meryl count k=15 output asm.meryl asm.fasta
	meryl print greater-than distinct=0.9998 asm.meryl >  asm_repetitive_k15.meryl
	winnowmap -W asm_repetitive_k15.meryl reference_genome.fasta asm.fasta > 	asm_vs_reference.paf

Gap-closing

We run a gap-closing tool YaGCloser with the following command:

	# map the assembly to the reads and generate a sorted bam file
	minimap2 --secondary=yes -ax map-pb asm_chrlevel.fasta reads.fasta | samtools view -hSb - | samtools sort - > reference.sorted.bam

	# index the bam file
	samtools index reference.sorted.bam

	# get gaps coordinate of the assembly
	detgaps asm_chrlevel.fasta > asm_chrlevel.gaps.bed

	# identify potential gaps to be filled
	python ./yagcloser/yagcloser.py -g asm_chrlevel.fasta reads.fasta -a reference.sorted.bam -b asm_chrlevel.gaps.bed -o yagcloser_output -f 100 -mins 10 -s  asm_chrlevel

Mitogenome assembly

Mitochondrial DNA is circular and multiple copies exist in the cell, it is often assembled with multiple duplication in one long contig, hence we utilised MitoHifi, a program that assembles the mitochondrial genome (mtDNA) using a reference sequence. To assemble the mtDNA, we used the following command:

	singularity exec --bind /usr/lib/locale/ docker://ghcr.io/marcelauliano/mitohifi:master 	mitohifi.py -c asm.fasta -f OP921772.1.fasta -g OP921772.1.gb -t 8 -o 2

Assembly assessment

To assess the completeness score of the assembly, we used BUSCO with the following parameters:

	busco -c 8 -i asm.fasta -l ./mammalia_odb10 -o busco_swamp_buffalo -m genome –offline

Alignment percent identity

For identifying percentage identity of the mitochondrial sequence, we use Blast+ with the following parameters:

	blastn -db ref.fa -query query.fa -outfmt 7 -evalue 1e-10 -perc_identity 90 -out 	blast_output.blnm7

Repeat Analysis

To analyse the repeats of the genome assembly, we used RepeatMasker with the default command:

	RepeatMasker -species "Bubalus bubalis" -dir output asm.fasta

Telomere and centromere

To determine the number of the telomeric repeat (TTAGGG)n in the assemblies, we used tidk with the following command:

	tidk search -s TTAGGG -o output --dir output_dir -e tsv asm.fasta
	tidk plot -o output_plot -t repeat_windows.tsv

For the centromeric repeats, we determined the satellite/ctr repeat families identified by RepeatMasker and analysed it using a custom Rscript. We then extract the span of the centromeric region using samtools to find tandem repeats using Tandem Repeat Finder (trf). The following parameters were used to find the monomer:

	trf centromere.fa 2 7 7 80 10 50 2000 -d -m -l 10
Next, we looked for the number of copies of the tandem repeats using HiCAT with the following argument:

	python HiCAT.py -i centromere.fa -t tandem_repeat.fa -o outfile –threads 2

Phylogeny and estimation of divergence

We used the salmonid_pipeline to get the orthogroups and single-copy genes of the species. After concatenating the single-copy genes, a phylogenetic tree was made using the following command with IQtree. This parameter also gives divergence estimation using LSD2 method.

iqtree -s bovine_8sp.txt -nt AUTO -B 1000 --date datefile_8sp.txt --date-tip 0 --date-ci 100

To estimate divergence time using PAML (mcmctree), we used the following configuration file:

	#mcmctree config file

	seed = -1
	seqfile = ../../data/bovine_8sp.txt
	treefile = ../../data/bovine_8sp.treefile
	mcmcfile = mcmc.txt
	outfile = out.txt

	ndata = 1
	seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
	usedata = 2    * 0: no data (prior); 1:exact likelihood;
	                      * 2:approximate likelihood; 3:out.BV (in.BV)
	clock = 2    * 1: global clock; 2: independent rates; 3: correlated rates
	RootAge = '<1.0'  * safe constraint on root age, used if no fossil for root.

	model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
	alpha = 0.5  * alpha for gamma rates at sites
	ncatG = 5    * No. categories in discrete gamma

	cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

	BDparas = 1 1 0   * birth, death, sampling
	kappa_gamma = 6 2     * gamma prior for kappa
	alpha_gamma = 1 1     * gamma prior for alpha

	rgene_gamma = 2 40 1   * gammaDir prior for rate for genes
	sigma2_gamma = 1 10 1   * gammaDir prior for sigma^2     (for clock=2 or 3)

	print = 1   * 0: no mcmc sample; 1: everything except branch rates 2: everything
	burnin = 20000
	sampfreq = 200
	nsample = 20000


Structural variations and SNPs

To determine the structural variations and SNPs with the water buffalo genome assemblies, we first removed the gaps using seqtk:

	seqtk cutN -n 3 ref.fa > ref_no_gaps.fa
	seqtk cutN -n 3 qry.fa > qry_no_gaps.fa

We then aligned the genomes to the reference genome using nucmer alignment with the following parameters:

	nucmer --maxmatch -l 100 -c 500  ref_no_gaps.fa qry_no_gaps.fa -p out -t 24

Next, we use the delta file to identify structural variants greater than 50 bp but no less than 10,000 bp with Assemblytics:

	Assemblytics out.delta out_assemblytics 10000 50 10000

The bedfile output were then analysed using a custom Rscript to remove duplicated calls and assess the length and number of structural variants.
