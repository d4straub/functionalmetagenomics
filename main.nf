//DEFAULT PARAMS

// mandatory
params.input = false // --input "data/*_R{1,2}.fastq.gz"

// singularity exec 'https://depot.galaxyproject.org/singularity/metaphlan:4.0.6--pyhca03a8a_0' metaphlan --install --index mpa_vOct22_CHOCOPhlAnSGB_202212 --bowtie2db metaphlan_database
params.metaphlan_database = false

// wget http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_DEMO_diamond_v201901b.tar.gz
// mkdir nuc_database
// tar -xzvf DEMO_chocophlan.v201901_v31.tar.gz -C nuc_database
// --humann_nuc_database "humann_dbs/nuc_database/*"
params.humann_nuc_database = false

// wget http://huttenhower.sph.harvard.edu/humann_data/chocophlan/DEMO_chocophlan.v201901_v31.tar.gz
// mkdir prot_database
// tar -xzvf uniref90_DEMO_diamond_v201901b.tar.gz -C prot_database
params.humann_prot_database = false

// optional DBs
params.kneaddata_genome_db = false
params.kneaddata_transcriptome_db = false
params.kneaddata_rrna_db = false

// optional
params.regroup_options = "--groups uniref90_rxn" // " "
params.regroup_customdb = false // "map_ko_uniref90.txt.gz"
params.rename_options = "--names metacyc-rxn" // "--names kegg-orthology"

// skip or stop at steps
params.skip_fastqc = false
params.stop_kneaddata = false
params.stop_metaphlan = false
params.stop_humann = false

//CHANNELS

if(!params.input){exit 1, log.info "Cannot find any input files. Example: --input \"data/*.fasta\""}

//PROCESSES

process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    container 'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Make list of old name and new name pairs to use for renaming in the bash while loop
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${prefix}.${reads.extension}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${prefix}_${index + 1}.${entry.extension}" ] }
    def rename_to = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ old_name, new_name -> new_name }.join(' ')
    """
    printf "%s %s\\n" $rename_to | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done

    fastqc \\
        $args \\
        --threads $task.cpus \\
        $renamed_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
}

process MULTIQC {
    label 'process_single'

    container 'https://depot.galaxyproject.org/singularity/multiqc:1.17--pyhdfd78af_0'

    input:
    path  multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    """
    multiqc \\
        --force \\
        $args \\
        $config \\
        $extra_config \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}

process KNEADDATA {
    tag "$meta.id"
    label 'process_high'

    container 'https://depot.galaxyproject.org/singularity/kneaddata:0.12.0--pyhdfd78af_1'

    input:
    tuple val(meta), path(reads) // only paired-end allowed
    path(genome_db)
    path(transcriptome_db)
    path(rrna_db)

    output:
    tuple val(meta), path("output/${meta.id}.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")                     , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}"
    def genome = genome_db ? '--reference-db ' + genome_db : ''
    def transcriptome = transcriptome_db ? '--reference-db ' + transcriptome_db : ''
    def rrna = rrna_db ? '--reference-db ' + rrna_db : ''
    """
    # this is needed with '--trimmomatic' because the 'trimmomatic.jar' is otherwise not found
    cp /usr/local/share/trimmomatic/trimmomatic.jar ./

    kneaddata \\
        --trimmomatic ./ \\
        --threads $task.cpus \\
        --input1 ${reads[0]} \\
        --input2 ${reads[1]} \\
        --output output \\
        --cat-final-output \\
        --output-prefix ${prefix} \\
        --log ${prefix}.log \\
        $genome \\
        $transcriptome \\
        $rrna \\
        $args
    gzip -n output/*.fastq
    """
}

process KNEADDATA_SUMMARY {
    label 'process_medium'

    container 'https://depot.galaxyproject.org/singularity/kneaddata:0.12.0--pyhdfd78af_1'

    input:
    path(input, stageAs: 'kneaddata_logs/*')

    output:
    path("kneaddata_read_count_table.tsv") , emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    kneaddata_read_count_table \\
        --input kneaddata_logs \\
        --output kneaddata_read_count_table.tsv
    """
}

process METAPHLAN_METAPHLAN {
    tag "$meta.id"
    label 'process_high'

    container 'https://depot.galaxyproject.org/singularity/metaphlan:4.0.6--pyhca03a8a_0'

    input:
    tuple val(meta), path(input)
    path metaphlan_db_latest

    output:
    tuple val(meta), path("*_profile.txt")   ,                emit: profile
    tuple val(meta), path("*.biom")          ,                emit: biom
    tuple val(meta), path('*.bowtie2out.txt'), optional:true, emit: bt2out
    path "versions.yml"                      ,                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_type = "$input" =~ /.*\.(fastq|fq)/ ? "--input_type fastq" : "$input" =~ /.*\.(fasta|fna|fa)/ ? "--input_type fasta" : "$input".endsWith(".bowtie2out.txt") ? "--input_type bowtie2out" : "--input_type sam"
    def input_data  = ("$input_type".contains("fastq")) && !meta.single_end ? "${input[0]},${input[1]}" : "$input"
    def bowtie2_out = "$input_type" == "--input_type bowtie2out" || "$input_type" == "--input_type sam" ? '' : "--bowtie2out ${prefix}.bowtie2out.txt"
    """
    BT2_DB=`find -L "${metaphlan_db_latest}" -name "*rev.1.bt2*" -exec dirname {} \\;`
    BT2_DB_INDEX=`find -L ${metaphlan_db_latest} -name "*.rev.1.bt2*" | sed 's/\\.rev.1.bt2.*\$//' | sed 's/.*\\///'`

    metaphlan \\
        --nproc $task.cpus \\
        $input_type \\
        $input_data \\
        $args \\
        $bowtie2_out \\
        --bowtie2db \$BT2_DB \\
        --index \$BT2_DB_INDEX \\
        --biom ${prefix}.biom \\
        --output_file ${prefix}_profile.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """
}

process METAPHLAN_MERGEMETAPHLANTABLES {
    label 'process_single'

    container 'https://depot.galaxyproject.org/singularity/metaphlan:4.0.6--pyhca03a8a_0'

    input:
    path(profiles)

    output:
    path("merged_profiles.txt") , emit: txt
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    merge_metaphlan_tables.py \\
        $args \\
        -o merged_profiles.txt \\
        ${profiles}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """
}

process HUMANN_HUMANN {
    tag "$meta.id"
    label 'process_high'

    container 'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0'

    input:
    tuple val(meta), path(reads), path(profile)
    path(nuc_database, stageAs: 'nuc_database/*')
    path(prot_database, stageAs: 'prot_database/*')

    output:
    tuple val(meta), path("*_genefamilies.tsv.gz") , emit: genefamilies
    tuple val(meta), path("*_pathabundance.tsv.gz"), emit: pathabundance
    tuple val(meta), path("*_pathcoverage.tsv.gz") , emit: pathcoverage
    tuple val(meta), path("*.log")                 , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    humann \\
        --threads $task.cpus \\
        --input $reads \\
        $args \\
        --nucleotide-database nuc_database \\
        --protein-database prot_database \\
        --taxonomic-profile $profile \\
        --output out \\
        --output-basename $prefix \\
        --o-log ${prefix}.log \\
        --remove-temp-output
    gzip -n out/*.tsv
    cp out/${prefix}_genefamilies.tsv.gz ${prefix}_genefamilies.tsv.gz
    cp out/${prefix}_pathabundance.tsv.gz ${prefix}_pathabundance.tsv.gz
    cp out/${prefix}_pathcoverage.tsv.gz ${prefix}_pathcoverage.tsv.gz
    """
}

process HUMANN_JOIN {
    label 'process_low'

    container 'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0'

    input:
    path('input/*')

    output:
    path("humann_genefamilies.tsv.gz"), emit: genefamilies
    path("humann_pathcoverage.tsv.gz"), emit: pathcoverage
    path("humann_pathabundance.tsv.gz"), emit: pathabundance
    path("*.log"), emit: log

    script:
    """
    humann_join_tables \\
        --input input \\
        --output humann_genefamilies.tsv \\
        --file_name genefamilies \\
        > join_genefamilies.log

    humann_join_tables \\
        --input input \\
        --output humann_pathcoverage.tsv \\
        --file_name pathcoverage \\
        > join_pathcoverage.log

    humann_join_tables \\
        --input input \\
        --output humann_pathabundance.tsv \\
        --file_name pathabundance \\
        > join_pathabundance.log

    gzip -n *.tsv
    """
}

process HUMANN_REGROUP {
    tag "$meta.id"
    label 'process_medium'

    container 'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0'

    input:
    tuple val(meta), path(input)
    val options
    path(customdb)

    output:
    tuple val(meta), path("*_regroup.tsv.gz"), emit: regroup
    tuple val(meta), path("*.log"), emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}"
    def cmd_customdb = customdb ? "--custom $customdb" : ""
    """
    humann_regroup_table \\
        --input $input \\
        --output ${prefix}_regroup.tsv \\
        $cmd_customdb \\
        $options \\
        > ${prefix}_regroup.log
    gzip -n ${prefix}_regroup.tsv
    """
}

process HUMANN_RENORM {
    tag "$meta.id"
    label 'process_medium'

    container 'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0'

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_renorm.tsv.gz"), emit: renorm
    tuple val(meta), path("*.log"), emit: log

    script:
    def prefix = "${meta.id}"
    """
    humann_renorm_table \\
        --input $input \\
        --output ${prefix}_renorm.tsv \\
        --units cpm \\
        --update-snames \\
        >${prefix}_renorm.log
    gzip -n ${prefix}_renorm.tsv
    """
}

process HUMANN_RENAME {
    tag "$meta.id"
    label 'process_medium'

    container 'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0'

    input:
    tuple val(meta), path(input)
    val options

    output:
    tuple val(meta), path("*_rename.tsv.gz"), emit: rename
    tuple val(meta), path("*.log"), emit: log

    script:
    def prefix = "${meta.id}"
    """
    humann_rename_table \\
        --input $input \\
        --output ${prefix}_rename.tsv \\
        $options \\
        > ${prefix}_rename.log

    gzip -n ${prefix}_rename.tsv
    """
}

process HUMANN_SPLITSTARTIFIEDTABLE {
    tag "$meta.id"
    label 'process_medium'

    container 'https://depot.galaxyproject.org/singularity/humann:3.8--pyh7cba7a3_0'

    input:
    tuple val(meta), path(input, stageAs: 'input/*')

    output:
    tuple val(meta), path("*_unstratified.tsv.gz"), emit: unstratified
    tuple val(meta), path("*_stratified.tsv.gz"), emit: stratified
    tuple val(meta), path("*.log"), emit: log

    script:
    def prefix = "${meta.id}"
    """
    humann_split_stratified_table \\
        --input $input \\
        --output . \\
        > ${input.baseName}_stratified.log

    gzip -n *.tsv
    """
}

//PIPELINE

workflow {
    // get files, paired-end required
    // after kneaddata (which expects here paired-end data, the data will be collated into one file, hence "single_end")
    Channel
        .fromFilePairs( params.input, size: 2 )
        .ifEmpty { error("Cannot find path file ${params.input}") }
        .map { name, reads ->
                def meta = [:]
                meta.id           = name
                meta.single_end   = true
                [ meta, reads ] }
        .set { ch_reads }

    // get DBs
    kneaddata_genome_db = params.kneaddata_genome_db ? file(params.kneaddata_genome_db) : []
    kneaddata_transcriptome_db = params.kneaddata_transcriptome_db ? file(params.kneaddata_transcriptome_db) : []
    kneaddata_rrna_db = params.kneaddata_rrna_db ? file(params.kneaddata_rrna_db) : []
    humann_nuc_database = params.humann_nuc_database ? file(params.humann_nuc_database) : []
    humann_prot_database = params.humann_prot_database ? file(params.humann_prot_database) : []
    metaphlan_database = params.metaphlan_database ? file(params.metaphlan_database) : []
    regroup_customdb = params.regroup_customdb ? file(params.regroup_customdb) : []

    // Process

    // 1 QC visualisations
    FASTQC ( params.skip_fastqc ? Channel.empty() : ch_reads )
    MULTIQC ( params.skip_fastqc ? Channel.empty() : FASTQC.out.zip.collect{it[1]}.ifEmpty([]), [], [], [] )

    // 2 Preprocessing
    KNEADDATA ( params.stop_kneaddata ? Channel.empty() : ch_reads, kneaddata_genome_db, kneaddata_transcriptome_db, kneaddata_rrna_db )
    KNEADDATA.out.log
        .map { meta, files -> files }
        .collect()
        .set { ch_kneaddata_summary }
    KNEADDATA_SUMMARY ( ch_kneaddata_summary.dump(tag: 'ch_kneaddata_summary') )
    ch_reads_preprocessed = KNEADDATA.out.reads

    // 3 Phylogenetic profiles
    METAPHLAN_METAPHLAN ( params.stop_metaphlan ? Channel.empty() : ch_reads_preprocessed, metaphlan_database )
    METAPHLAN_METAPHLAN.out.profile
        .map { meta, files -> files }
        .collect()
        .set { ch_metaphlan_summary }
    METAPHLAN_MERGEMETAPHLANTABLES ( ch_metaphlan_summary )

    // 4 Functional profiles
    ch_reads_preprocessed
        .join( METAPHLAN_METAPHLAN.out.profile )
        .set { ch_humann }
    HUMANN_HUMANN ( params.stop_humann ? Channel.empty() : ch_humann, humann_nuc_database, humann_prot_database )
    HUMANN_HUMANN.out.genefamilies
        .mix( HUMANN_HUMANN.out.pathabundance)
        .mix( HUMANN_HUMANN.out.pathcoverage)
        .map { meta, files -> [ files ] }
        .collect()
        .set { ch_join }
    HUMANN_JOIN ( ch_join )

    // 5 Postprocessing
    HUMANN_JOIN.out.genefamilies
        .map { read ->
                def meta = [:]
                meta.id = "genefamilies"
                [ meta, read ] }
        .set{ch_humann_genefamilies}
    HUMANN_JOIN.out.pathabundance
        .map { read ->
            def meta = [:]
            meta.id = "pathabundance"
            [ meta, read ] }
        .set{ch_humann_pathabundance}
    HUMANN_JOIN.out.pathcoverage
        .map { read ->
            def meta = [:]
            meta.id = "pathcoverage"
            [ meta, read ] }
        .set{ch_humann_pathcoverage}

    // this is only for genefamilies and pathabundance, not for pathcoverage
    HUMANN_RENORM ( ch_humann_genefamilies.mix(ch_humann_pathabundance) )
    HUMANN_RENORM.out.renorm
        .branch {
            // it[0] is meta map!
            genefamilies: it[0].id == "genefamilies"
            pathabundance: it[0].id == "pathabundance"
        }
        .set { ch_renormed }
    ch_renormed.genefamilies.set { ch_renormed_genefamilies }
    ch_renormed.pathabundance.set { ch_renormed_pathabundance }

    // this is only for genefamilies
    HUMANN_REGROUP ( ch_renormed_genefamilies, params.regroup_options, regroup_customdb )
    HUMANN_RENAME ( HUMANN_REGROUP.out.regroup, params.rename_options )

    // this is for all types
    ch_to_unstratify = HUMANN_RENAME.out.rename.mix(ch_renormed_pathabundance).mix(ch_humann_pathcoverage)
    HUMANN_SPLITSTARTIFIEDTABLE ( ch_to_unstratify )
}
