#!/usr/bin/env nextflow

//Clean main to start coding

import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
nextflow.enable.dsl = 2
nextflow.preview.recursion=true 

params.inputfile = "runMT-human_ont.mmi"


// seq_data_folder = "/home/charlotte/heackaton/HEK_HeLa_Experiment_folder/HEK_HeLa/20230323_1311_MN32609_FAQ51879_03abc106/fastq_pass/"
// output_dir = "~/out_dir"
// metadata = null
// barcoded = null
// genome_gtf = "/home/charlotte/heackaton/gencode.v43.primary_assembly.annotation.chr20.gtf"
// bed_file = null
// genome_fasta = null
// genome_index = "runMT-human_ont.mmi"
// transcriptome_fasta = null

include { convert_gtf_to_df } from './lib/initialize.nf'

process Minimap {
    label "wftemplate"

       publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )
    input:
    tuple val(ID), path(fastq)
    file(index)

    output: 
    path "${ID}.out${task.index}.bam" 
 
  
    script:

    """
    minimap2 --MD -ax splice -uf -k14 -t ${params.threads} ${params.inputfile} ${params.fastq} | samtools view -hbS -F 3844 | samtools sort > ${ID}.out${task.index}.bam
    """

}
process FeatureCount {
    label "wftemplate"
                   

       publishDir(
        path: "${params.outdir}",
        mode: 'copy',
    )
    input:
    path bam

    output: 
    path "${bam.simpleName}.merged_fc${task.index}.csv", emit: csv
 
  
    script:

    """
    samtools index ${bam}
	featureCounts -a "${params.gtf}" -F 'GTF' -L -T ${params.threads} -o ${bam.simpleName}.merged_fc${task.index}.csv ${bam}
    """
    
}

process progressive_stats {
   

    input: 
        path output
    output:
       //path("all_stats.${task.index}")
       path "merged_fc${task.index}.csv"
    script:
        def new_input = output instanceof BlankSeparatedList ? output.first() : output
       def state = output instanceof BlankSeparatedList ? output.last() : "/home/charlotte/heackaton/Charlotte/~/out_dir/test.csv"
       def output = "feature_counts_latest_${task.index}.csv"
    """
    Rscript ${projectDir}/script.r "${new_input}" "${state}"
    mv /home/charlotte/heackaton/Charlotte/~/out_dir/feature_counts_latest.csv  merged_fc${task.index}.csv
    """
}



///home/charlotte/heackaton/HEK_HeLa_Experiment_folder/HEK_HeLa/20230323_1311_MN32609_FAQ51879_03abc106/fastq_pass
//initialise workflow and hand over input
workflow {
    convert_gtf_to_df(file("${params.genome_gtf}"))
	data = Channel      
    .watchPath("${params.seq_data_folder}/*/*.fastq.gz", 'create,modify')      
    .until { file->file.name == 'STOP.fastq.gz' }
        .map { tuple(it.parent.name, it ) }
        .set { samples }
         samples.view()


    //test_output	= Minimap(samples, params.inputfile)
    test_output	= Minimap(samples, file(params.inputfile))
    // Here you have the output channel as a collection
    test_output.view()

    output	= FeatureCount(test_output)

   // def merged = FeatureCount.out.csv.collectFile(keepHeader: true).view()
    // Here you have the output channel as a collection
    test_output.view()

    progressive_stats.scan(output)
  
}
//messages to display once the workflow has completed

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'Finish Correctly :)' : 'uh oh sth went wrong' }"
}

//Rscript script.r "${new_input}" "${state}"