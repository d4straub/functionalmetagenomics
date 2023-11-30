# functionalmetagenomics

Assembly-free functional metagenomics pipeline with KneadData and HUMAnN3

## Software requirement

- Java
- nextflow
- singularity / apptainer

## Usage

### Input data

Input data has to be Illumina paired-end data in compressed fastq files (e.g. `*.fq.gz`).
Download data into a folder, here proposed `data`.

With `--input 'data/*_R{1,2}_L001.fastq.gz'` any paired files that follow the regex will be used (if just one file is present, it will be ignored). Any string matching `*` will be used as sample name, e.g. the paired files `data/A_B_C_R1_L001.fastq.gz` & `data/A_B_C_R2_L001.fastq.gz` will get the sample name `A_B_C`.
If the data is in nested folders, such as `data/*/<files>`, use `--input 'data/*/*_R{1,2}_L001.fastq.gz'` (the folder name will be ignored).

### Prepare databases

#### KneadData

These are databases to remove human genome, transcriptome or rRNA sequences. All of those are optional.
Refer to the [KneadData documentation](https://huttenhower.sph.harvard.edu/kneaddata/).

```bash
singularity pull kneaddata_0-12-0.sif 'https://depot.galaxyproject.org/singularity/kneaddata:0.12.0--pyhdfd78af_1'
singularity exec kneaddata_0-12-0.sif kneaddata_database --download human_genome bowtie2 human_genome
singularity exec kneaddata_0-12-0.sif kneaddata_database --download human_transcriptome bowtie2 human_transcriptome
singularity exec kneaddata_0-12-0.sif kneaddata_database --download ribosomal_RNA bowtie2 ribosomal_RNA
```

#### MetaPhlAn4

This is a database for taxonomic classification, it is required.
Refer to the [MetaPhlAn documentation](https://huttenhower.sph.harvard.edu/metaphlan/).

```bash
singularity pull metaphlan_4-0-6.sif  'https://depot.galaxyproject.org/singularity/metaphlan:4.0.6--pyhca03a8a_0'
singularity exec metaphlan_4-0-6.sif metaphlan --install --index mpa_vOct22_CHOCOPhlAnSGB_202212 --bowtie2db metaphlan_database
```

#### HUMAnN3

These are databases for functional classification. Other databases are available (e.g. uniref50_diamond might be used).
The DB `utility_mapping` is required when regrouping to KO or other identifiers.
Refer to the [HUMAnN3 documentation](https://huttenhower.sph.harvard.edu/humann/).

```bash
singularity pull humann_3-8.sif  'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0'
singularity exec humann_3-8.sif humann_databases --download chocophlan full humann_dbs --update-config no
singularity exec humann_3-8.sif humann_databases --download uniref uniref90_diamond humann_dbs --update-config no
singularity exec humann_3-8.sif humann_databases --download utility_mapping full humann_dbs --update-config no
```

### Run the pipeline

```bash
NXF_VER=23.10.0 nextflow run d4straub/functionalmetagenomics -r 1.0.0 \
  --input 'data/*_R{1,2}_L001.fastq.gz' \
  --kneaddata_genome_db 'human_genome' \
  --kneaddata_transcriptome_db 'human_transcriptome' \
  --kneaddata_rrna_db 'ribosomal_RNA' \
  --metaphlan_database 'metaphlan_database' \
  --humann_nuc_database 'humann_dbs/chocophlan/*' \
  --humann_prot_database 'humann_dbs/uniref/*' \
  --outdir result
```

All `--kneaddata*` parameters are optional but recommended.

> **Warning:**
> Use a config to specify your cluster specifications, e.g. the here available `cfc.config` by appending to the above command `-c cfc.config`.

If re-grouping to KO is desired, append to the upper command:

```bash
  --regroup_options " " \
  --regroup_customdb "humann_dbs/utility_mapping/map_ko_uniref90.txt.gz" \
  --rename_options "--names kegg-orthology"
```

Additional parameters are
- `--skip_fastqc` -> skip QC = FastQC & MultiQC
- `--stop_kneaddata` -> stop before KneadData step (i.e. only QC = FastQC & MultiQC)
- `--stop_metaphlan` -> stop before MetaPhlAn step (i.e. QC and preprocessing)
- `--stop_humann` --> stop before HUMAnN (i.e. only QC, preprocessing, and taxonomic community profiling)

## Output

This section describes only the most important output files.

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.

[KneadData](https://huttenhower.sph.harvard.edu/kneaddata/) is a tool designed to perform quality control on metagenomic sequencing data.

- `kneaddata/`
  - `kneaddata_read_count_table.tsv`:tab-separated file with read counts during read filtering for all samples
 
[MetaPhlAn](https://huttenhower.sph.harvard.edu/metaphlan/) is a computational tool for profiling the composition of microbial communities.

- `metaphlan/`
  - `merged_profiles.txt`: tab-separated file with taxa and abundances for all samples

[HUMAnN3](https://huttenhower.sph.harvard.edu/humann/) is a method for efficiently and accurately profiling the abundance of microbial metabolic pathways and other molecular functions. For more details on what each file contains, refer to the manual of HUMAnN3.

- `humann/`
  - `0_joined_raw`: tables incuding all samples, direct output of HUMAnN
    - `humann_genefamilies.tsv.gz`: functional annotations and their abundances
    - `humann_pathabundance.tsv.gz`: pathway abundances
    - `humann_pathcoverage.tsv.gz`: pathway coverage
  - `4_final`: normalized (CPM), summarized, and named tables incuding all samples, stratified (per taxa) or unstratified (taxa independent)
    - `genefamilies_rename_{stratified/unstratified}.tsv.gz`: functional annotations and their abundances
    - `pathabundance_renorm_{stratified/unstratified}.tsv.gz`: pathway abundances
    - `humann_pathcoverage_{stratified/unstratified}.tsv.gz`: pathway coverage

## Credits

This pipeline was originally written by Daniel Straub ([@d4straub](https://github.com/d4straub)) for use at the [Quantitative Biology Center (QBiC)](http://www.qbic.life).

Code was inspired by nf-core/ampliseq (doi: 10.5281/zenodo.1493841) ([Straub et al., 2020](https://doi.org/10.3389/fmicb.2020.550420)) of the nf-core collection of workflows ([Ewels et al., 2020](https://dx.doi.org/10.1038/s41587-020-0439-x)).

## Citations

Please cite all tools, such as

- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)
> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)
  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)
  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

- [KneadData](https://huttenhower.sph.harvard.edu/kneaddata/)

- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
  > Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/faq.shtml)
  > Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.

- [MetaPhlAn](https://huttenhower.sph.harvard.edu/metaphlan/)
  > Aitor Blanco-Miguez, Francesco Beghini, Fabio Cumbo, Lauren J. McIver, Kelsey N. Thompson, Moreno Zolfo, Paolo Manghi, Leonard Dubois, Kun D. Huang, Andrew Maltez Thomas, Gianmarco Piccinno, Elisa Piperni, Michal Punčochář, Mireia Valles-Colomer, Adrian Tett, Francesca Giordano, Richard Davies, Jonathan Wolf, Sarah E. Berry, Tim D. Spector, Eric A. Franzosa, Edoardo Pasolli, Francesco Asnicar, Curtis Huttenhower, Nicola Segata. Extending and improving metagenomic taxonomic profiling with uncharacterized species with MetaPhlAn 4. bioRxiv 2022.08.22.504593; doi: https://doi.org/10.1101/2022.08.22.504593 

- [HUMAnN3](https://huttenhower.sph.harvard.edu/humann/)
  > Beghini F, McIver LJ, Blanco-Míguez A, Dubois L, Asnicar F, Maharjan S, Mailyan A, Manghi P, Scholz M, Thomas AM, Valles-Colomer M, Weingart G, Zhang Y, Zolfo M, Huttenhower C, Franzosa EA, Segata N. Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. Elife. 2021 May 4;10:e65088. doi: 10.7554/eLife.65088. PMID: 33944776; PMCID: PMC8096432.
