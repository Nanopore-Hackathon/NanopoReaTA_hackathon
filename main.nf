#!/usr/bin/env nextflow

//Clean main to start coding

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/ingress' //This is how you can integrate subworkflows


OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE") //Empty file to be used as gapfiller if some files can not be provided due to issues

process x {
    label "wftemplate"
    publishDir ( //Use publish dir to write outputs to Folders in your filesystem
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
    
    output:

    script:

    """
    #
    """
}

// workflow module
workflow pipeline {
    take:

    main:

}

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
