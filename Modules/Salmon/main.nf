process salmon_annotation {
    memory "4 GB"
    input:
        val string from salmon_channel
        val string_array from salmon_channel2

    output:
        val string into salmon_annotation_merge
        val string_array into salmon_annotation_merge2

    script:
        if (string != "") { 
            data_string = ""
            for (i in 0..params.sample_names.size()-1){
                data_string = data_string + params.run_dir + params.sample_names.get(i) + "/" + "bam_files_transcripts" + "/" + "full" + "/" +  " "
            }
            salmon_annotation_running.value= 1 
        
        bam_string_transcripts=""
        for(i in 0..string_array.size()-1) {
        if (string_array.size() > 0) {
            element = string_array.get(i)
            if (params.barcoded == 1){
            basis = element.getParent().getName()
            }
            else {
            basis = element.getParent().getParent().getParent().getName()
            }
            filename=element.getName()
            converted_filename =  filename.minus(".${params.suffix}")
            outname=params.run_dir + basis + "/bam_files_transcripts/full/" + converted_filename + ".bam" + " "
            bam_string_transcripts=bam_string_transcripts + outname
            }
            }
        }

        if (string != "")
            """
            function samtoolsParallel2 {
                if [ ! \$(ls \${1}| wc -l) -eq 0 ] 
                then
                    sample=\$1
                    run_dir=\$2
                    par_basis=\$(basename \$(dirname \$(dirname \${sample})))
                    if [ ! -f \${run_dir}bam_transcriptome_merged/\${par_basis}.bam ]
                    then
                        samtools merge \${run_dir}\${par_basis}/salmon/all.bam \${sample}*.bam -f --threads ${params.threads} -c -p || echo "\${sample}" >> \${run_dir}/error_logs/merge_transcriptome_all_samtools_error.log
                        cp \${run_dir}\${par_basis}/salmon/all.bam \${run_dir}bam_transcriptome_merged/\${par_basis}.bam 
                        samtools index \${run_dir}bam_transcriptome_merged/\${par_basis}.bam || echo "Indexing failed" >> \${run_dir}/error_logs/indexing_transcriptome_all_error.log
                    else
                    bam_files_to_merge=\${run_dir}\${par_basis}/salmon/all.bam" "
                    for i in ${bam_string_transcripts}
                    do
                    if [[ "\$(basename \$(dirname \$(dirname \$(dirname \${i}))))" == "\${par_basis}" ]]
                        then
                        bam_files_to_merge=\${bam_files_to_merge}\${i}" "
                        fi
                    done
                    samtools merge \${run_dir}bam_transcriptome_merged/\${par_basis}.bam \${bam_files_to_merge} -f --threads ${params.threads} -c -p || echo "\${bam_files_to_merge}" >> \${run_dir}/error_logs/merge_transcriptome_few_error.log
                    cp \${run_dir}bam_transcriptome_merged/\${par_basis}.bam \${run_dir}\${par_basis}/salmon/all.bam
                    samtools index \${run_dir}bam_transcriptome_merged/\${par_basis}.bam --threads ${params.threads} || echo "Indexing failed" >> \${run_dir}/error_logs/indexing_transcriptome_few_error.log
                    fi
                fi
            }
            export -f samtoolsParallel2 
            start=\$(date +%s)
            parallel -v -u --env samtoolsParallel2 --no-notice -j 5 samtoolsParallel2 ::: ${data_string} ::: ${params.run_dir}
            end=\$(date +%s)
            time=\$(echo "\$((\$end-\$start))")
            echo "samtools_transcriptome,${iteration.value},\$time" >> ${params.run_dir}processing_time_table.csv
            start=\$(date +%s)
            for i in ${data_string}
            do 
                basis=\$(basename \$(dirname \$(dirname \${i})))
                if [ ! \$(ls \${i}| wc -l) -eq 0 ]
                then
                salmon quant -t ${params.transcriptome_fasta} -l A -a ${params.run_dir}\${basis}/salmon/all.bam -o ${params.run_dir}\${basis}/salmon/ -g ${params.genome_gtf} -p ${params.threads} || echo "\$i " >> ${params.run_dir}error_logs/salmon_annotation_fail.log
                fi
            done
            end=\$(date +%s)
            time=\$(echo "\$((\$end-\$start))")
            echo "salmon,${iteration.value},\$time" >> ${params.run_dir}processing_time_table.csv 
            """
        else
            """
            echo ""
            """
}

process merge_salmon_annotation {
    memory "4 GB"
    input:
        val string from salmon_annotation_merge
        val string_array from salmon_annotation_merge2
    
    output:
        val string into salmon_annotation_add_names
        val string_array into salmon_annotation_add_names2

    script:
        println "Salmon merging all files"  
        data_string = ""
        for (i in 0..params.sample_names.size()-1){
        data_string = data_string + params.run_dir + params.sample_names.get(i) + "/" + " "
        }
        if (string != "") 
            """
            python ${params.script_dir}merge_all_salmon.py ${data_string} ${params.run_dir}salmon_merged_tpm.csv ${params.run_dir}salmon_merged_absolute.csv || echo "${string}" >> ${params.run_dir}error_logs/salmon_merge_error.log
            """
        else
            """
            echo ""
            """
}

process salmon_annotation_done{
    memory "4 GB"
    input:
        val string from salmon_annotation_add_names
        val string_array from salmon_annotation_add_names2

    output:
        val done into salmon_done_channel
    
    exec:
        if (string != ""){
            salmon_annotation_running.value = 0
        }
        done = 1
}