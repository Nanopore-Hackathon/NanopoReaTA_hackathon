import groovy.json.JsonBuilder
@Grab('com.xlson.groovycsv:groovycsv:1.1')
import static com.xlson.groovycsv.CsvParser.parseCsv
import java.io.File;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                 FeatureCount Annotation                                                                 //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process count_features{
    memory "4 GB"
    input:
    val bam_string
    val string_array 
    output:
    val fc_string
    val string_array

    script:
    fc_string = ""
    fc_string_splice = ""
    data_string = ""
    for(i in 0..params.sample_names.size()-1){
       if(string_array.size() > 0){
           basis = params.sample_names.get(i)
           outname= params.output_dir + basis + "/merged_fc/merged_fc.csv" + " "
           fc_string= fc_string + outname
           data_string = data_string + params.sample_names.get(i) + " "
       }
    } 
    if (bam_string != "")
    """
    start=\$(date +%s)
    for basis in ${data_string}
    do
        featureCounts -a "${params.genome_gtf}" -F 'GTF' -L -T ${params.threads} -o ${params.output_dir}\${basis}/merged_fc/merged_fc.csv ${params.output_dir}\${basis}/all_gene.bam || echo \${basis} >> ${params.output_dir}error_logs/featureCounts_error.log
    done
    end=\$(date +%s)
    time=\$(echo "\$((\$end-\$start))")
    echo "featureCounts,${iteration.value},\$time" >> ${params.output_dir}processing_time_table.csv
    """
    else
    """
    echo ""
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                 Salmon Annotation Transcriptome                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process salmon_annotation{
    memory "4 GB"
    input:
        val string 
        val string_array 

    output:
        val string 
        val string_array 

    script:
    if (string != ""){ 
        data_string = ""
        for (i in 0..params.sample_names.size()-1){
	        data_string = data_string + params.output_dir + params.sample_names.get(i) + "/" + "bam_files_transcripts" + "/" + "full" + "/" +  " "
        }
        salmon_annotation_running.value= 1 
    
    bam_string_transcripts=""
    for(i in 0..string_array.size()-1){
       if (string_array.size() > 0){
           element = string_array.get(i)
           if (params.barcoded == 1){
           basis = element.getParent().getName()
           }
           else {
               basis = element.getParent().getParent().getParent().getName()
           }
           filename=element.getName()
           converted_filename =  filename.minus(".${params.suffix}")
           outname=params.output_dir + basis + "/bam_files_transcripts/full/" + converted_filename + ".bam" + " "
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
        parallel -v -u --env samtoolsParallel2 --no-notice -j 5 samtoolsParallel2 ::: ${data_string} ::: ${params.output_dir}
        end=\$(date +%s)
        time=\$(echo "\$((\$end-\$start))")
        echo "samtools_transcriptome,${iteration.value},\$time" >> ${params.output_dir}processing_time_table.csv
        start=\$(date +%s)
        for i in ${data_string}
        do 
            basis=\$(basename \$(dirname \$(dirname \${i})))
            if [ ! \$(ls \${i}| wc -l) -eq 0 ]
            then
            salmon quant -t ${params.transcriptome_fasta} -l A -a ${params.output_dir}\${basis}/salmon/all.bam -o ${params.output_dir}\${basis}/salmon/ -g ${params.genome_gtf} -p ${params.threads} || echo "\$i " >> ${params.output_dir}error_logs/salmon_annotation_fail.log
            fi
        done
        end=\$(date +%s)
        time=\$(echo "\$((\$end-\$start))")
        echo "salmon,${iteration.value},\$time" >> ${params.output_dir}processing_time_table.csv 
        """
    else
        """
        echo ""
        """
}
