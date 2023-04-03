# Setup and initial analyses

Create a project directory from which run commands, for instance:

```sh
mkdir Super-sinks
cd Super-sinks
```

Download reads from BioProject **PRJNA803387** in `../Reads_PF_TC/raw_reads`. Trim reads with Trim Galore v0.6.6, indicating with the --basename flag the strain name as in the BioSample record (PF_* and TC_* for plasmid free and transconjugant strains, respectively):

```sh
trim_galore --quality 20 --illumina --length 50 --fastqc --basename WTCHG_<strain_name> --output_dir ../Reads_PF_TC/trimmed/ --paired ../Reads_PF_TC/raw_reads/<fastq1> ../Reads_PF_TC/raw_reads/<fastq2>
```

From BioProject **PRJNA626430**, download the raw reads and genome assemblies of the 200 *E. coli* and *Klebsiella* spp. strains carrying pOXA-48 listed in DelaFuente et al. 2022 (https://doi.org/10.1038/s41559-022-01908-7), Supplementary Data 1. Trim reads with Trim Galore v0.6.4 as above and store trimmed reads in `../Reads_WT_pOXA-48/trimmed`. Assemblies should be stored in `../Assemblies_WT_pOXA-48`.

Assembly of PF and TC strains with SPAdes v3.15.2 and assembly quality assessment with QUAST v5.0.2:

```sh
for fq1 in ../Reads_PF_TC/trimmed/*val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	strain=${fq1::-12}
	strain=${strain:29}
	spades.py --isolate --cov-cutoff auto -o ../Assemblies_PF_TC/$strain -1 $fq1 -2 $fq2
	cp ../Assemblies_PF_TC/$strain/contigs.fasta ../Assemblies_PF_TC/$strain"_contigs.fasta"
done

quast.py -o ../Assemblies_PF_TC/quast/ --glimmer ../Assemblies_PF_TC/*fasta
```

TC were confirmed to correspond to their PF strain by checking mutations and plasmid content with Snippy v4.6.0 and ABRicate v1.0.1 (PlasmidFinder database), respectively. Snippy was also used to confirm isogenicity of the acquired pOXA-48 plasmids, using the sequence of pOXA-48_K8 (MT441554), located in `../Closed_sequences/plasmids/`,  as reference. Multilocus sequence typing was performed with MLST v2.21.0.

```sh
for file in ../Assemblies_PF_TC/PF*contigs.fasta
do
	strain=${file:23}
	strain=${strain::-14}
	fq1="../Reads_PF_TC/trimmed/WTCHG_TC_"$strain"_val_1.fq.gz"
	fq2="../Reads_PF_TC/trimmed/WTCHG_TC_"$strain"_val_2.fq.gz"
	# checking correspondance PF-TC
	snippy --report --outdir ../Assemblies_PF_TC/analysis/variants_snippy/ref_PF_map_TC_$strain --ref $file --R1 $fq1 --R2 $fq2
	# checking isogenicity with pOXA-48_K8
	snippy --report --outdir ../Assemblies_PF_TC/analysis/variants_snippy/ref_pOXA-48_K8_map_TC_$strain --ref ../Closed_sequences/plasmids/pOXA-48_K8.fasta --R1 $fq1 --R2 $fq2
done

mkdir ../Assemblies_PF_TC/analysis/abricate/
for file in ../Assemblies_PF_TC/*contigs.fasta
do
	name=${file:20}
	name=${name::-14}
	abricate --db plasmidfinder $file  ../Assemblies_PF_TC/analysis/abricate/$name
done
abricate --summary ../Assemblies_PF_TC/analysis/abricate/* > ../Assemblies_PF_TC/analysis/summary_plasmids.tsv

mlst ../Assemblies_PF_TC/*contigs.fasta | cut -f1,3 > ../Assemblies_PF_TC/analysis/ST.tsv
```


Mash distance phylogenies were constructed for the 25 E. coli and 25 Klebsiella spp. PF strains from their whole-genome assemblies using mashtree v1.2.0 with a bootstrap of 100:

```sh
mkdir ../Assemblies_PF_TC/phylogeny_mash_PF
mashtree_bootstrap.pl --reps 100 ../Assemblies_PF_TC/PF_EC*_contigs.fasta > ../Assemblies_PF_TC/phylogeny_mash_PF/mashtree_btrp_Ec.dnd
mashtree_bootstrap.pl --reps 100 ../Assemblies_PF_TC/PF_K*_contigs.fasta > ../Assemblies_PF_TC/phylogeny_mash_PF/mashtree_btrp_Kleb.dnd
```

# Understanding differences in acquired resistance levels

## Plasmid copy number

```sh
mkdir plasmid_copy_number
cd plasmid_copy_number
```

Download the **mapping_script.sh**, **coverage2pcn.sh** and **strain_species.tsv** files in the current directory. Execute `./mapping_script.sh` to map the trimmed reads to their respective genome assemblies using BWA MEM v0.7.17.

Plasmid replicons are identified in the assemblies using ABRicate v1.0.1 with the PlasmidFinder database. Next, the results are filtered to obtain the contig containing the pOXA-48 replicon.

```sh
mkdir plasmids_abricate
# pOXA-48 carriers
for contigs in ../../Assemblies_WT_pOXA-48/*.fasta
do
	strain=${contigs:28}
	strain=${strain::-14}
	abricate --db plasmidfinder $contigs > ./plasmids_abricate/$strain".tsv"
done
# Transconjugants
for contigs in ../../Assemblies_PF_TC/TC*.fasta
do
	strain=${contigs:23}
	strain=${strain::-14}
	abricate --db plasmidfinder $contigs > ./plasmids_abricate/$strain".tsv"
done
# Filtering results
for file in ./plasmids_abricate/*.tsv
do
	strain=${file:20}
	strain=${strain::-4}
	cat $file | grep "IncL/M(pOXA-48)" > ./plasmids_abricate/$strain"_filt.tsv"
done
```

Execute `./coverage2pcn.sh` to calculate the pOXA-48 copy number (PCN), estimated as the ratio of pOXA-48/chromosome median coverage. In this study, the median coverage of pOXA-48 is calculated from the contig containing the IncL replicon, and the median coverage of the chromosome is calculated from the first three largest contigs (NODE_1, NODE_2 and NODE_3), confirmed to be chromosomal through BLASTn to the NCBI database (see **blastn_chr123\*.txt** files). Additionally, the script also calculates the mean pOXA-48 and chromosomal coverage, as well as median and mean chromosomal coverage from only the first, largest chromosomal contig, and the PCN ratios for all these calculations. Before running the script, be sure to install GNU Datamash v1.4 and to modify the script as indicated in its code. You can redirect the output to a tsv file with `./coverage2pcn.sh > PCN_WT_pOXA-48.tsv` or `./coverage2pcn.sh > PCN_TC.tsv`. The R script **PCN_plots.R** contains the code used to generate the figures.

Files **PF_TC_MICs.csv** (MICs of tested antibiotics of PF and TC strains) and **TC_MICs_FC_PCN.csv** (fold change in MIC of TC strains and PCN) are provided with the **stats_PCN_MIC.R** script, which computes statistics of Figures 1A and S3.

## Antimicrobial resistance genes

The genomes of the PF strains were searched for antimicrobial resistance genes using ABRicate v1.0.1 with the ResFinder database.

```sh
cd ..
mkdir amr_mic_fold_change
for file in ../Assemblies_PF_TC/PF*_contigs.fasta
do
	strain=${file:20}
	strain=${strain::-14}
	abricate --db resfinder $file > amr_mic_fold_change/$strain".tsv"
done
abricate --summary amr_mic_fold_change/* > amr_mic_fold_change/summary_amr.tsv
```

The **summary_amr.tsv** file was edited to replace dots for 0s (absence of the gene). The edited file is provided along with the **AMR_plots.R** script to generate figures.


# Understanding differences in pOXA-48 acquisition by conjugation

Copy the fasta files of the 20 recipients and 3 donors to `./contigs_recipients_donors` and annotate genomes with Prokka v1.14.6:

```sh
mkdir contigs_recipients_donors
mkdir analysis_recipients
# Recipients
cp ../Assemblies_PF_TC/PF_EC02_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_EC03_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_EC07_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_EC10_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_EC12_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_EC16_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_EC19_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_EC21_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_EC22_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_EC25_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_KPN01_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_KPN02_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_KPN09_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_KPN10_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_KPN11_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_KPN12_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_KPN14_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_KPN16_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_KQ01_contigs.fasta contigs_recipients_donors/
cp ../Assemblies_PF_TC/PF_KQ04_contigs.fasta contigs_recipients_donors/
# Donors
cp ../Assemblies_WT_pOXA-48/C165_contigs.fasta contigs_recipients_donors/donor_C165_contigs.fasta
cp ../Assemblies_WT_pOXA-48/K93_contigs.fasta contigs_recipients_donors/donor_K93_contigs.fasta
# Annotation
for file in contigs_recipients_donors/*contigs.fasta
do
	strain=${file:26}
	strain=${strain::-14}
	prokka -o contigs_recipients_donors/prokka/$strain --prefix $strain $file
	cp contigs_recipients_donors/prokka/$strain/$strain".faa" contigs_recipients_donors/
	cp contigs_recipients_donors/prokka/$strain/$strain".gff" contigs_recipients_donors/
done
# For donor E. coli β3914 download the protein fasta sequences of E. coli K-12 (NC_000913.3) in contigs_recipients_donors/donor_MG1655.faa
```

### Plasmid incompatibility

Plasmid replicons were detected with ABRicate v1.0.1 with the PlasmidFinder database:

```sh
mkdir analysis_recipients/plasmids_abricate
for file in contigs_recipients_donors/PF*contigs.fasta
do
	strain=${file:26}
	strain=${strain::-14}
	abricate --db plasmidfinder $file > analysis_recipients/plasmids_abricate/$strain".tsv"
done
abricate --summary analysis_recipients/plasmids_abricate/* > analysis_recipients/plasmids_abricate/summary_reps.tsv
```

The **summary_reps.tsv** file was edited to replace dots for 0s (absence of the gene). The edited file is provided along with the **recipients_plots.R** script to generate figures.


### Restriction-modification (RM) systems

Download the **REBASE gold_standard database** of protein sequences from http://rebase.neb.com/rebase/rebase.seqs.html in the next directory:

```sh
mkdir analysis_recipients/defense_REBASE
```

The protein sequences of the recipients and donors were blasted (BLASTp v2.9.0) against the REBASE database:

```sh
makeblastdb -in analysis_recipients/defense_REBASE/protein_gold_seqs.txt -dbtype prot
for file in contigs_recipients_donors/*faa
do
	strain=${file:26}
	strain=${strain::-4}
	blastp -query $file -db protein_gold_seqs.txt > analysis_recipients/defense_REBASE/$strain".tsv"
done
```

The BLASTp results were parsed with the **parser_blastp.py** script, which incorporates enzyme information from the REBASE database and filters matches by percentage of identity and alignment based on the criteria from Oliveira et al. 2016 (https://doi.org/10.1073/pnas.1603257113). Briefly, hits with type I and IIG RM systems were kept if the percentage of identity and alignment was 100%, since it was observed that recognition sites varied between enzymes with only a few mismatches. For type II systems, a threshold of 70% identity and 90% alignment was used, as it was previously reported that these systems with similarity over 55% generally share the same target motifs. For type III and type IV this threshold was higher, as reported in the same study, and was set at 90% identity and 90% alignment. After filtering, some proteins had hits with more than one enzyme. For these cases, the best hit or the hit with the enzyme of the corresponding organism (*E. coli* or *Klebsiella* spp.) was kept.

RM systems were also detected in the proteomes of the recipients and donors using the DefenseFinder webserver (accessed February 2022, results stored in `analysis_recipients/defense_finder_web/<strain>`). Then, the REBASE hits and the RM systems detected by DefenseFinder were unified by protein ID. Enzymes or RM systems present in all strains were discarded, since they would not explain differences in conjugation frequencies between isolates. Finally, only complete RM systems that were not present in the donor strains were retained for each recipient-donor combination, since the donor would be conferring protection to pOXA-48 against equivalent systems of the recipient. Complete type I systems comprise a restriction enzyme (REase), a methylase (MTase) and a specificity (S) subunit. Type II systems include at least a REase and a MTase. No type III systems were detected, and type IV systems are normally composed of one or two proteins. Donor and recipient systems were regarded as similar if they shared the same recognition sequence or if a BLASTp alignment of the enzymes followed the previous criteria of identity and alignment percentage per type of system. The same method was applied for comparing the final systems between recipient strains, and similar systems within each type were given the same letter code identifying the RM system subtype. The final presence/absence tables (**RM_sys.csv** and **RM_sys_sum.csv**) are provided along with the **recipients_plots.R** script to generate figures.

### CRISPR arrays

The nucleotide sequences of the 20 recipients were submited to the DefenseFinder webserver to detect CRISPR arrays. Results were downloaded to `./analysis_recipients/defense_finder_web_CRISPR/<strain>`. The nucleotide sequences of the CRISPR arrays are located in `./analysis_recipients/defense_finder_web_CRISPR/<strain>/crisprcas-out/rawCRISPRs.fna`, which contain spacers and direct repeats (dr). The spacers within these sequences are annotated in `./analysis_recipients/defense_finder_web_CRISPR/<strain>/crisprcas-out/GFF/<contig_file.gff>`. The script **parser_crisprs.py** reads the GFF file of a contig containing CRISPR arrays (provided manually after inspecting the rawCRISPRs.fna file) and extracts the spacer sequences. Therefore, it can be run as following, concatenating the spacers in the spacers.fasta file per strain.

```sh
./analysis_recipients/defense_finder_web_CRISPR/parser_crisprs.py ./analysis_recipients/defense_finder_web_CRISPR/<strain>/crisprcas-out/GFF/<contig_file.gff> >> ./analysis_recipients/defense_finder_web_CRISPR/<strain>/spacers.fasta
```

Next, the spacer sequences were blasted (BLASTn v2.9.0) against the sequence of pOXA-48_K8 (MT441554), located in `../Closed_sequences/plasmids/`. Results are stored in blastn_spacers_pOXA-48_K8.tsv in each strain's directory.

```sh
for file in analysis_recipients/defense_finder_web_CRISPR/*/spacers.fasta
do
	dir=${file::-13}
	echo -e "qseqid\tsseqid\tpident\tlength\tqlen\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > $dir"blastn_spacers_pOXA-48_K8.tsv"
	blastn -query $file -db ../Closed_sequences/plasmids/pOXA-48_K8.fasta -outfmt "6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore" -task blastn-short >> $dir"blastn_spacers_pOXA-48_K8.tsv"
done
```

The script **parser_blastn.py** calculates the percentage of alignment of each hit and filters by percentage of alignment (tunable inside script). The best alignment was 63% in strain PF_KPN02, therefore CRISPR systems are not targeting pOXA-48.

```sh
for dir in analysis_recipients/defense_finder_web_CRISPR/*; do analysis_recipients/defense_finder_web_CRISPR/parser_blastn.py $dir; done
```

### Other defense systems

Other defense systems were identified from the previous DefenseFinder results (stored in `analysis_recipients/defense_finder_web/<strain>`) and from running PADLOC v1.1.0 (database v1.3.0):

```sh
for file in contigs_recipients_donors/PF*faa
do
	name=${file:26}
	name=${name::-4}
	padloc --faa $file --gff $name".gff" --outdir analysis_recipients/defense_padloc/
done
```

The DefenseFinder and PADLOC results (csv file) were merged manually by protein ID, discarding Cas systems. The final presence/absence table (**defensesys.csv**) is provided along with the **recipients_plots.R** script to generate figures.

### Type VI secretion systems

Secretion systems were searched using the TXSScan models of MacSyFinder v1.0.5. Only complete systems were considered.

```sh
for faa in contigs_recipients_donors/PF*faa
do
	name=${fasta:26}
	name=${fasta::-4}
	macsyfinder T6SSi T6SSii T6SSiii --db-type ordered_replicon --replicon-topology linear -d ~/Software/macsy_models/TXSS/definitions/ -p ~/Software/macsy_models/TXSS/profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_recipients/t6ss/$name
done
```

### Capsules

*Klebsiella* capsules were identified for all recipient strains with Kaptive v2.0.0. The presence of *Klebsiella* capsules was assigned when the confidence level of the matches was above "good".

```sh
for fasta in contigs_recipients_donors/PF*contigs.fasta
do
	kaptive.py -a $fasta -k ~/Software/Kaptive/reference_database/Klebsiella_k_locus_primary_reference.gbk -o analysis_recipients/capsules_kaptive_Kpn_Ec/ -t 1
done
```

Other capsule systems were searched using MacSyFinder v1.0.5. Only systems reported as complete were considered in following analyses.

```sh
mkdir analysis_recipients/capsules_capsulefinder/
mkdir analysis_recipients/capsules_capsulefinder/wzy
mkdir analysis_recipients/capsules_capsulefinder/abc
mkdir analysis_recipients/capsules_capsulefinder/groupiv_e
mkdir analysis_recipients/capsules_capsulefinder/groupiv_f
mkdir analysis_recipients/capsules_capsulefinder/groupiv_s
mkdir analysis_recipients/capsules_capsulefinder/pga
mkdir analysis_recipients/capsules_capsulefinder/syn_cps3
mkdir analysis_recipients/capsules_capsulefinder/syn_has

for faa in contigs_recipients_donors/PF*faa
do
	name=${faa:26}
	name=${name::-4}
	macsyfinder Wzy_stricte --db-type ordered_replicon --replicon-topology linear -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_recipients/capsules_capsulefinder/wzy/$name
	macsyfinder ABC --db-type ordered_replicon --replicon-topology linear -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_recipients/capsules_capsulefinder/abc/$name
	macsyfinder GroupIV_e_stricte --db-type ordered_replicon --replicon-topology linear -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_recipients/capsules_capsulefinder/groupiv_e/$name
	macsyfinder GroupIV_f --db-type ordered_replicon --replicon-topology linear -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_recipients/capsules_capsulefinder/groupiv_f/$name
	macsyfinder GroupIV_s_stricte --db-type ordered_replicon --replicon-topology linear -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_recipients/capsules_capsulefinder/groupiv_s/$name
	macsyfinder PGA --db-type ordered_replicon --replicon-topology linear -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_recipients/capsules_capsulefinder/pga/$name
	macsyfinder Syn_cps3 --db-type ordered_replicon --replicon-topology linear -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_recipients/capsules_capsulefinder/syn_cps3/$name
	macsyfinder Syn_has --db-type ordered_replicon --replicon-topology linear -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_recipients/capsules_capsulefinder/syn_has/$name
done
```

To assess the similarity of the different capsular system sets (*Klebsiella* capsules and/or other capsular systems encoded by a strain) between strains, a weighted gene repertoire relatedness (**wGRR**) score was calculated as in de Sousa et al. 2020 (https://doi.org/10.1038/s41396-020-0726-z), using a python script developed by Rocha's lab. The pipeline requires as input a single multifasta file including the protein sequences of the capsular operons identified for all strains (**caps_sets_prots.faa**). For *Klebsiella* capsules, each operon nucleotide sequence was annotated with Prokka v1.14.6 (default parameters) to obtain the protein sequences. For other capsular systems, the protein sequences within the borders of the systems reported by CapsuleFinder were concatenated in a single fasta. The script output (**wGRR_genomes_Twoway.txt**) is provided along with the **heatmap_wGRR.R** script to generate the wGRR heatmap.

## Analysis of *Klebsiella* capsules in a database

The NCBI RefSeq database of high quality complete non-redundant *K. pneumoniae* (n=730) and *E. coli* (n=1,585) genomes was retrieved from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/ on March 2021 and downloaded in `database_refseq/`:

```sh
mkdir database_refseq
mkdir analysis_databases
```

Strain names were assigned as in Table S3, where the first four letters indicate species (ESCO for *E. coli* and KLPN for *K. pneumoniae*). Several fasta files were downloaded from each strain. The following is an example of files downloaded for one strain, ESCO001.0321.00010:
* ESCO001.0321.00010.C001.fst --> complete nucleotide sequence of the chromosome
* ESCO001.0321.00010.P002.fst, ESCO001.0321.00010.P003.fst --> complete nucleotide sequences of each plasmid
* ESCO001.0321.00010.prt --> all protein sequences (chromosome and plasmids)

### Identifying pOXA-48 carrying strains

pOXA-48-carrying strains were identified in the RefSeq database using three approaches. (i) Plasmid replicons and (ii) antimicrobial resistance genes were identified with ABRicate v1.0.1. The names of the strains having a pOXA-48 replicon and/or the *bla*<sub>OXA-48</sub> gene were stored in different files.

```sh
mkdir -p analysis_databases/plasmids_resistance_refseq/abricate_plas/
mkdir analysis_databases/plasmids_resistance_refseq/abricate_res/

for fasta in database_refseq/*fst
do
	name=${fasta:16}
	name=${name::-4}
	abricate --db plasmidfinder $fasta > analysis_databases/plasmids_resistance_refseq/abricate_plas/$name".tsv"
	abricate --db resfinder $fasta > analysis_databases/plasmids_resistance_refseq/abricate_res/$name".tsv"
done
# pOXA-48 replicon
grep "pOXA-48" analysis_databases/plasmids_resistance_refseq/abricate_plas/* | cut -f 2 | sort | uniq > analysis_databases/plasmids_resistance_refseq/abricate_names_pOXA-48.txt
# blaOXA-48 gene
grep "blaOXA-48" analysis_databases/plasmids_resistance_refseq/abricate_res/* | cut -f 2 | sort | uniq > analysis_databases/plasmids_resistance_refseq/abricate_names_blaOXA-48.txt
```

(iii) The sequence of pOXA-48_K8 (MT441554) was blasted (BLASTn v2.9.0) against the RefSeq database. Hits with >95% identity and >10000 bp alignment with pOXA-48_K8 were filtered and the names of the strains passing the filteres were stored.

```sh
for fasta in database_refseq/*fst
do
	name=${fasta:16}
	name=${name::-4}
	makeblastdb -in $fasta -dbtype nucl
	blastn -query ../Closed_sequences/plasmids/pOXA-48_K8.fasta -db $fasta -outfmt 6 > analysis_databases/plasmids_resistance_refseq/blastn/$name".tsv"
done
# Filtering
awk '{if ($3 > 95 && $4 > 10000) { print } }' analysis_databases/plasmids_resistance_refseq/blastn/* | cut -f 2 | sort | uniq > analysis_databases/plasmids_resistance_refseq/blastn_names_pOXA-48.txt
```

Then, pOXA-48-carriage was defined as (i) having the pOXA-48 replicon, (ii) having the *bla*<sub>OXA-48</sub> gene and (iii) passing the filters of identity and alignment to pOXA-48_K8:

```sh
comm -12 analysis_databases/plasmids_resistance_refseq/blastn_names_pOXA-48.txt analysis_databases/plasmids_resistance_refseq/abricate_names_pOXA-48.txt | comm -12 - analysis_databases/plasmids_resistance_refseq/abricate_names_blaOXA-48.txt > analysis_databases/plasmids_resistance_refseq/common_names_pOXA-48_plasmid.txt
cat analysis_databases/plasmids_resistance_refseq/common_names_pOXA-48_plasmid.txt | sed 's/.....$//' > analysis_databases/plasmids_resistance_refseq/common_names_pOXA-48_strain.txt
```

Thus, the file **common_names_pOXA-48_strain.txt** contains the names of the strains that carry pOXA-48. There are only 52 strains out of 2315 that carry pOXA-48, of which 6 are *E. coli* and 46 *K. pneumoniae*.

### Plasmids in the database

Next lines calculate the number of strains that don't carry plasmids and that carry increasing number of plasmids. The maximum number of different plasmids carried by one strain is 13. At 6 plasmids, the number of strains falls below 100. Thus, for creating Fig. 4C, analyses were performed with increasing number of plasmids until ≥6 plasmids.

```sh
# Number of strains carring X number of plasmids
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 1) {print}}' | awk -F " " '{print $2}' | wc -l # 0 plasmids
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 2) {print}}' | awk -F " " '{print $2}' | wc -l
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 3) {print}}' | awk -F " " '{print $2}' | wc -l
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 4) {print}}' | awk -F " " '{print $2}' | wc -l
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 5) {print}}' | awk -F " " '{print $2}' | wc -l
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 6) {print}}' | awk -F " " '{print $2}' | wc -l
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 7) {print}}' | awk -F " " '{print $2}' | wc -l
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 8) {print}}' | awk -F " " '{print $2}' | wc -l
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 9) {print}}' | awk -F " " '{print $2}' | wc -l
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 10) {print}}' | awk -F " " '{print $2}' | wc -l
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 11) {print}}' | awk -F " " '{print $2}' | wc -l
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 12) {print}}' | awk -F " " '{print $2}' | wc -l
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 13) {print}}' | awk -F " " '{print $2}' | wc -l
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 14) {print}}' | awk -F " " '{print $2}' | wc -l # 13 plasmids
# Saving strain names with increasing number of plasmids
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 >= 2) {print}}' | awk -F " " '{print $2}' > database_refseq/names_strains_with_plasmids.txt
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 1) {print}}' | awk -F " " '{print $2}' > database_refseq/names_strains_with_eq0_plasmids.txt
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 2) {print}}' | awk -F " " '{print $2}' > database_refseq/names_strains_with_eq1_plasmids.txt
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 3) {print}}' | awk -F " " '{print $2}' > database_refseq/names_strains_with_eq2_plasmids.txt
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 4) {print}}' | awk -F " " '{print $2}' > database_refseq/names_strains_with_eq3_plasmids.txt
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 5) {print}}' | awk -F " " '{print $2}' > database_refseq/names_strains_with_eq4_plasmids.txt
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 == 6) {print}}' | awk -F " " '{print $2}' > database_refseq/names_strains_with_eq5_plasmids.txt
ls database_refseq/*fst | grep -E ".C0|.P0" | sed 's/.........$//' | sed 's/database_refseq\///' | sort | uniq -c | awk -F" " '{if ($1 >= 7) {print}}' | awk -F " " '{print $2}' > database_refseq/names_strains_with_eqmore6_plasmids.txt
cat database_refseq/names_strains_with_plasmids.txt <(ls database_R-GNOSIS/*/*fasta | cut -d"/" -f 3 | cut -d"." -f 1) > analysis_databases/names_strains_with_plasmids_RefSeq+R-GNOSIS_all.txt
```

File **names_strains_with_plasmids.txt** contains names of strains with >=1 plasmids in RefSeq.

### Identifying capsules in the database

Kaptive v2.0.0 was run as before to identify *Klebsiella* capsules:

```sh
# RefSeq E. coli
for fasta in database_refseq/ESCO*fst
do
	kaptive.py -a $fasta -k ~/Software/Kaptive/reference_database/Klebsiella_k_locus_primary_reference.gbk -o analysis_databases/capsules_kaptive_refseq_Ec/ -t 1
done
# RefSeq K. pneumoniae
for fasta in database_refseq/KLPN*fst
do
	kaptive.py -a $fasta -k ~/Software/Kaptive/reference_database/Klebsiella_k_locus_primary_reference.gbk -o analysis_databases/capsules_kaptive_refseq_Kpn/ -t 1
done
```

CapsuleFinder (MacSyFinder v1.0.5) was used to identify other capsule systems to exclude non-capsulated bacteria from following analyses.

```sh
# RefSeq E. coli
mkdir analysis_databases/capsules_capsulefinder_refseq_Ec/
mkdir analysis_databases/capsules_capsulefinder_refseq_Ec/wzy
mkdir analysis_databases/capsules_capsulefinder_refseq_Ec/abc
mkdir analysis_databases/capsules_capsulefinder_refseq_Ec/groupiv_e
mkdir analysis_databases/capsules_capsulefinder_refseq_Ec/groupiv_f
mkdir analysis_databases/capsules_capsulefinder_refseq_Ec/groupiv_s
mkdir analysis_databases/capsules_capsulefinder_refseq_Ec/pga
mkdir analysis_databases/capsules_capsulefinder_refseq_Ec/syn_cps3
mkdir analysis_databases/capsules_capsulefinder_refseq_Ec/syn_has
for faa in database_refseq/ESCO*.*.[0-9][0-9][0-9][0-9][0-9].prt
do
	name=${faa:16}
	name=${name::-4}
	macsyfinder -w 3 Wzy_stricte --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Ec/wzy/$name
	macsyfinder -w 3 ABC --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Ec/abc/$name
	macsyfinder -w 3 GroupIV_e_stricte --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Ec/groupiv_e/$name
	macsyfinder -w 3 GroupIV_f --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Ec/groupiv_f/$name
	macsyfinder -w 3 GroupIV_s_stricte --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Ec/groupiv_s/$name
	macsyfinder -w 3 PGA --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Ec/pga/$name
	macsyfinder -w 3 Syn_cps3 --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Ec/syn_cps3/$name
	macsyfinder -w 3 Syn_has --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Ec/syn_has/$name
done

# RefSeq K. pneumoniae
mkdir analysis_databases/capsules_capsulefinder_refseq_Kpn/
mkdir analysis_databases/capsules_capsulefinder_refseq_Kpn/wzy
mkdir analysis_databases/capsules_capsulefinder_refseq_Kpn/abc
mkdir analysis_databases/capsules_capsulefinder_refseq_Kpn/groupiv_e
mkdir analysis_databases/capsules_capsulefinder_refseq_Kpn/groupiv_f
mkdir analysis_databases/capsules_capsulefinder_refseq_Kpn/groupiv_s
mkdir analysis_databases/capsules_capsulefinder_refseq_Kpn/pga
mkdir analysis_databases/capsules_capsulefinder_refseq_Kpn/syn_cps3
mkdir analysis_databases/capsules_capsulefinder_refseq_Kpn/syn_has
for faa in database_refseq/KLPN*.*.[0-9][0-9][0-9][0-9][0-9].prt
do
	name=${faa:16}
	name=${name::-4}
	macsyfinder -w 3 Wzy_stricte --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Kpn/wzy/$name
	macsyfinder -w 3 ABC --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Kpn/abc/$name
	macsyfinder -w 3 GroupIV_e_stricte --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Kpn/groupiv_e/$name
	macsyfinder -w 3 GroupIV_f --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Kpn/groupiv_f/$name
	macsyfinder -w 3 GroupIV_s_stricte --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Kpn/groupiv_s/$name
	macsyfinder -w 3 PGA --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Kpn/pga/$name
	macsyfinder -w 3 Syn_cps3 --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Kpn/syn_cps3/$name
	macsyfinder -w 3 Syn_has --db-type ordered_replicon -d ~/Software/macsy_models/capsulefinder/CapsuleFinder_models/Diderm_bacteria/ -p ~/Software/macsy_models/capsulefinder/CapsuleFinder_profiles/ --profile-suffix .hmm --sequence-db $faa -o analysis_databases/capsules_capsulefinder_refseq_Kpn/syn_has/$name
done
```

Strains were considered non-capsulated when CapsuleFinder did not find any complete system (the report file for each type of capsule was empty) and Kaptive results were below "good":

```sh
# RefSeq E. coli - all strains are capsulated
comm -23 <(find analysis_databases/capsules_capsulefinder_refseq_Ec/*/*/*report -size 0 | cut -d/ -f4 | sort | uniq -c | awk '{if ($1 > 7) print $2}') <(grep -E "Good|High|Very high|Perfect" analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | cut -f1 | cut -d_ -f 1 | sort | uniq | sed 's/.........$//')
# RefSeq K. pneumoniae - 15 strains are non-capsulated: KLPN001.0321.00003, KLPN001.0321.00259, KLPN001.0321.00260, KLPN001.0321.00261, KLPN001.0321.00335, KLPN001.0321.00389, KLPN001.0321.00390, KLPN001.0321.00391, KLPN001.0321.00394, KLPN001.0321.00440, KLPN001.0321.00470, KLPN001.0321.00520, KLPN001.0321.00599, KLPN001.0321.00725, KLPN001.0321.00727
comm -23 <(find analysis_databases/capsules_capsulefinder_refseq_Kpn/*/*/*report -size 0 | cut -d/ -f4 | sort | uniq -c | awk '{if ($1 > 7) print $2}') <(grep -E "Good|High|Very high|Perfect" analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | cut -f1 | cut -d_ -f 1 | sort | uniq | sed 's/.........$//')
```

### Testing associations

The Rmarkdown file **fisher_tests.Rmd** contains the code for building 2x2 contingency tables and running Fisher's exact tests. The association of enconding *Klebsiella* capsules *vs* other capsular systems and carrying pOXA-48 or plasmids is analyzed. Next lines calculate numbers necessary for building the contingency tables in R:

```sh
### Numbers for Fisher tests

# number of refseq Ec that have kleb caps
grep -E "Good|High|Very high|Perfect" analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | wc -l
# number of refseq Ec that do not have kleb caps
grep -E "None|Low" analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep ".C001" | wc -l
# number of refseq Ec that have kleb caps and pOXA-48
grep -f analysis_databases/plasmids_resistance_refseq/common_names_pOXA-48_strain.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
# number of refseq Ec that do not have kleb caps and have pOXA-48
grep -f analysis_databases/plasmids_resistance_refseq/common_names_pOXA-48_strain.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | wc -l
# number of refseq Ec that have kleb caps and plasmids
grep -f database_refseq/names_strains_with_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
# number of refseq Ec that do not have kleb caps and have plasmids
grep -f database_refseq/names_strains_with_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | wc -l

# number of refseq Kpn that have kleb caps
grep -E "Good|High|Very high|Perfect" analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | wc -l
# number of refseq Kpn that do not have kleb caps
grep -E "None|Low" analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
# number of refseq Kpn that have kleb caps and pOXA-48
grep -f analysis_databases/plasmids_resistance_refseq/common_names_pOXA-48_strain.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
# number of refseq Kpn that do not have kleb caps and have pOXA-48
grep -f analysis_databases/plasmids_resistance_refseq/common_names_pOXA-48_strain.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
# number of refseq Kpn that have kleb caps and plasmids
grep -f database_refseq/names_strains_with_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
# number of refseq Kpn that do not have kleb caps and have plasmids
grep -f database_refseq/names_strains_with_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l

### Numbers for graph odds ratio per number of plasmids (Fig. 4C)
## E. coli + K. pneumoniae
# have kleb caps and plasmids
grep -f database_refseq/names_strains_with_eq0_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq1_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq2_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq3_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq4_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq5_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eqmore6_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
# do not have kleb caps but have plasmids
grep -f database_refseq/names_strains_with_eq0_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
grep -f database_refseq/names_strains_with_eq1_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
grep -f database_refseq/names_strains_with_eq2_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
grep -f database_refseq/names_strains_with_eq3_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
grep -f database_refseq/names_strains_with_eq4_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
grep -f database_refseq/names_strains_with_eq5_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
grep -f database_refseq/names_strains_with_eqmore6_plasmids.txt analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
## E. coli
# have kleb caps and plasmids
grep -f database_refseq/names_strains_with_eq0_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq1_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq2_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq3_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq4_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq5_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eqmore6_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
# do not have kleb caps but have plasmids
grep -f database_refseq/names_strains_with_eq0_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | wc -l
grep -f database_refseq/names_strains_with_eq1_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | wc -l
grep -f database_refseq/names_strains_with_eq2_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | wc -l
grep -f database_refseq/names_strains_with_eq3_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | wc -l
grep -f database_refseq/names_strains_with_eq4_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | wc -l
grep -f database_refseq/names_strains_with_eq5_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | wc -l
grep -f database_refseq/names_strains_with_eqmore6_plasmids.txt analysis_databases/capsules_kaptive_refseq_Ec/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | wc -l
## K. pneumoniae
# have kleb caps and plasmids
grep -f database_refseq/names_strains_with_eq0_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq1_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq2_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq3_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq4_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eq5_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
grep -f database_refseq/names_strains_with_eqmore6_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "Good|High|Very high|Perfect" | wc -l
# do not have kleb caps but have plasmids
grep -f database_refseq/names_strains_with_eq0_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
grep -f database_refseq/names_strains_with_eq1_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
grep -f database_refseq/names_strains_with_eq2_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
grep -f database_refseq/names_strains_with_eq3_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
grep -f database_refseq/names_strains_with_eq4_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
grep -f database_refseq/names_strains_with_eq5_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l
grep -f database_refseq/names_strains_with_eqmore6_plasmids.txt analysis_databases/capsules_kaptive_refseq_Kpn/kaptive_results_table.txt | grep -E "None|Low" | grep ".C001" | grep -Ev "KLPN001.0321.00003|KLPN001.0321.00259|KLPN001.0321.00260|KLPN001.0321.00261|KLPN001.0321.00335|KLPN001.0321.00389|KLPN001.0321.00390|KLPN001.0321.00391|KLPN001.0321.00394|KLPN001.0321.00440|KLPN001.0321.00470|KLPN001.0321.00520|KLPN001.0321.00599|KLPN001.0321.00725|KLPN001.0321.00727" | wc -l

### Saving a necessary file for phylogenetic regression, which contains names of strains encoding Klebsiella capsules:
grep -E "Good|High|Very high|Perfect" analysis_databases/capsules_kaptive_refseq_*/kaptive_results_table.txt | cut -d":" -f2 | cut -f 1 | sed 's/.........$//g' > analysis_databases/names_strains_caps_RefSeq_all.txt
```

### Correcting for phylogenetic dependency

Significant interactions (P < 0.05, Fisher exact test) were then corrected for phylogenetic dependency, i.e. the tendency of closely related strains to share the same traits. For this, phylogenies were constructed from a set of conserved bacterial single copy genes. These are 139 gene markers present in >90% of bacterial species, described in Rinke et al. 2013 (https://doi.org/10.1038/nature12352, Suppl. Table 13). However, d’Humières et al. 2019 (https://doi.org/10.1038/s41598-019-47656-w) curated this list to remove genes present in phages, resulting in 128 HMM profiles specific to bacteria (Table S4). First, the Pfam-A HMM models (release 35) were downloaded and the 128 HMM profiles were extracted:

```sh
mkdir analysis_databases/hmm_gene_markers # Download here Pfam-A.hmm
cd analysis_databases/hmm_gene_markers # Analyses will be performed from this directory
# Extracting the 128 HMM profiles:
awk 'BEGIN{RS=ORS="//"}/PF03485|PF03484|PF01121|PF03772|PF06418|PF02768|PF00035|PF00889|PF01176|PF00113|PF03952|PF06574|PF03147|PF01687|PF02938|PF02527|PF00958|PF01025|PF01018|PF11987|PF04760|PF00707|PF05198|PF01715|PF06421|PF01795|PF02873|PF08529|PF01252|PF01195|PF00162|PF02912|PF03726|PF01416|PF02033|PF02132|PF00825|PF00687|PF00466|PF00298|PF03946|PF00572|PF00238|PF00252|PF01196|PF00828|PF00861|PF01245|PF00181|PF00453|PF00829|PF00237|PF00276|PF01016|PF00830|PF00831|PF03947|PF00297|PF01783|PF01632|PF00573|PF00281|PF00673|PF00347|PF03948|PF01281|PF00338|PF00411|PF00416|PF00312|PF00886|PF00366|PF01084|PF00203|PF00318|PF01649|PF00189|PF00163|PF00333|PF03719|PF01250|PF00177|PF00410|PF00380|PF00164|PF01782|PF01000|PF01193|PF04997|PF00623|PF04983|PF05000|PF04998|PF04563|PF04561|PF04565|PF10385|PF00562|PF04560|PF01765|PF02410|PF07499|PF01330|PF05491|PF03652|PF02773|PF02772|PF00584|PF03840|PF00344|PF02403|PF01668|PF02978|PF00763|PF02882|PF00121|PF03461|PF05698|PF05697|PF00750|PF01746|PF01509|PF02367|PF02130|PF12344|PF08459|PF10458|PF06071/' Pfam-A.hmm > bacterial_128_Pfam35.hmm
```

The protein fasta sequences of the 1,585 *E. coli* and 715 *K. pneumoniae* capsulated RefSeq strains were copied to `analysis_databases/hmm_gene_markers/refseq` (changing the extension from .prt to .faa). The protein fasta of the *Enterobacter cloacae* GCF_003204095.1 outgroup was also placed in `analysis_databases/hmm_gene_markers/refseq`, also with .faa extension. The HMM profiles were searched in the proteomes using HMMER v3.3 (option --cut_ga):

```sh
for faa in refseq/*faa
do
	name=${faa:7}
	name=${name::-4}
	hmmsearch --tblout refseq_hmm/$name".hmm128.txt" --cpu 20 --cut_ga bacterial_128_Pfam35.hmm $faa
done
```

Then, the top HMM hits were filtered by score using the cutoffs reported in Rinke et al. 2013 (https://doi.org/10.1038/nature12352, Suppl. Table 13), included in the **bacterial_128_Pfam35_list.tsv** file. This filtering step results in 127 gene markers present in >90% of strains:

```sh
# Get the top hit of each profile per proteome
for hit in refseq_hmm/*
do
    name=${hit::-4}
	awk '!x[$4]++' $hit | grep -v "#" > $name".top.txt"
done
# Filter by score cutoff
for hit in refseq_hmm/*top.txt
do
	name=${hit::-4}
	awk 'NR==FNR{a[$1]=$3;next} ($4 in a) && ($6 > a[$4]) {print $0}' bacterial_128_Pfam35_list.tsv $hit > $name".filt.txt"
done
# Number of strains (left column) having the corresponding number of filtered HMM hits (right column)
for hit in refseq_hmm/*top.filt.txt; do cat $hit | wc -l; done | sort | uniq -c
```

Protein sequences of each family (extracted and incorporated into a multifasta file with **hmmhits2fasta.py**) were aligned with MAFFT v7.453 (option --auto) and alignments were trimmed with trimAl v1.4.rev15. Then, alignments were concatenated with catfasta2phyml.pl (https://github.com/nylander/catfasta2phyml):

```sh
# Get the protein match of each HMM id
grep -E "ESCO|KLPN|WP_" refseq_hmm/*filt.txt | awk '{print $1"\t"$4}' | awk -F":" '{print $1"\t"$2}' | cut -f2,3 > prot_hmmid.tsv
# Split table by HMM id into different files
awk '{print >> $2".txt"}' prot_hmmid.tsv
# Remove the PF column
for pf in PF*; do sed -i 's/\tPF[0-9]*\.[0-9]*//g' $pf; done
# Get fasta sequences
./hmmhits2fasta.py

# MAFFT alignment and trimming with trimAl
mkdir refseq_alignments # will include Ec+Kpn+outgroup
mkdir refseq_alignments_Ec # will include Ec+outgroup
for faa in refseq_hmm/PF*faa
do
	name1=${faa::-4}
	name2=${name1:11}
	sed -n '/>ESCO/,/>KLPN/{/>KLPN/!p}' $faa > $name1"_Ec.faa"
	mafftout1="refseq_alignments/"$name2".mafft.aln"
	trimalout1="refseq_alignments/"$name2".mafft.trim.aln"
	mafftout2="refseq_alignments_Ec/"$name2".mafft.aln"
	trimalout2="refseq_alignments_Ec/"$name2".mafft.trim.aln"
	mafft --thread -1 --auto $faa > $mafftout1
	trimal -in $mafftout1 -out $trimalout1 -automated1
	mafft --thread -1 --auto $name1"_Ec.faa" > $mafftout2
	trimal -in $mafftout2 -out $trimalout2 -automated1
done

# Edit fasta headers to include only strain names
for file in refseq_alignments/*trim.aln
do
	name=${file::-4}
	awk '/^>/{print substr($1,1,19); next}{print $0}' $file | sed 's/WP_[0-9]*.[0-9]*/GCF_003204095.1/g' > $name".edit.aln"
done
for file in refseq_alignments_Ec/*trim.aln
do
	name=${file::-4}
	awk '/^>/{print substr($1,1,19); next}{print $0}' $file | sed 's/WP_[0-9]*.[0-9]*/GCF_003204095.1/g' > $name".edit.aln"
done
# Concatenate fastas
catfasta2phyml.pl -f -c refseq_alignments/*edit.aln > refseq_alignments/concat_all.mafft.trim.msa
catfasta2phyml.pl -f -c refseq_alignments_Ec/*edit.aln > refseq_alignments_Ec/concat_Ec.mafft.trim.msa
```

IQ-TREE v1.6.12 was used to infer two phylogenetic trees from the concatenated alignments (the first including all capsulated strains, and the second including only *E. coli* strains, both with the *E. cloacae* outgroup) with best evolutionary model selection and 1000 ultrafast bootstrap:

```sh
# 1st tree: Alignment has 2301 sequences with 47077 columns, 6248 distinct patterns
iqtree -s refseq_alignments/concat_all.mafft.trim.msa -nt AUTO -bb 1000 -m MFP
# 2nd tree: Alignment has 1586 sequences with 46491 columns, 4645 distinct patterns
iqtree -s refseq_alignments_Ec/concat_Ec.mafft.trim.msa -nt AUTO -bb 1000 -m MFP

cd ../.. # Back to the Super-sinks directory
mkdir analysis_databases/phylogenetic_regression_phyloglm
```

The phylogenetic trees (provided with the extension **.treefile**) are loaded to the R script **phylogenetic_regression.R** located in `analysis_databases/phylogenetic_regression_phyloglm/`. Trees are rooted and rescaled to a total hight of 1. The ordered tree nodes are exported to `analysis_databases/phylogenetic_regression_phyloglm/` under the names of treenodes_all.txt and treenodes_Ec.txt. These files, along with files generated before that contain the strain names that have the dependent (presence of pOXA-48) and independent (presence of *Klebsiella* capsule) traits, are fed to the **get_traits_table.py** script to obtain the traits tables:

```sh
# E. coli and K. pneumoniae from refseq carrying pOXA-48
analysis_databases/phylogenetic_regression_phyloglm/get_traits_table.py -n analysis_databases/phylogenetic_regression_phyloglm/treenodes_all.txt -d analysis_databases/plasmids_resistance_refseq/common_names_pOXA-48_strain.txt -i analysis_databases/names_strains_caps_RefSeq_all.txt > analysis_databases/phylogenetic_regression_phyloglm/traits_caps_pOXA-48_refseq_all.tsv
# E. coli and K. pneumoniae from refseq carrying plasmids
analysis_databases/phylogenetic_regression_phyloglm/get_traits_table.py -n analysis_databases/phylogenetic_regression_phyloglm/treenodes_all.txt -d database_refseq/names_strains_with_plasmids.txt -i analysis_databases/names_strains_caps_RefSeq_all.txt > analysis_databases/phylogenetic_regression_phyloglm/traits_caps_plasmids_refseq_all.tsv
# E. coli from RefSeq carrying plasmids
analysis_databases/phylogenetic_regression_phyloglm/get_traits_table.py -n analysis_databases/phylogenetic_regression_phyloglm/treenodes_Ec.txt -d database_refseq/names_strains_with_plasmids.txt -i analysis_databases/names_strains_caps_RefSeq_all.txt > analysis_databases/phylogenetic_regression_phyloglm/traits_caps_plasmids_refseq_Ec.tsv
```

These traits tables are read into the **phylogenetic_regression.R** script to perform phylogenetic logistic regression with the *phyloglm* v2.6.4 function, fitting a logistic MPLE model with 100 independent bootstrap replicates.
