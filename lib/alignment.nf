import groovy.json.JsonBuilder
@Grab('com.xlson.groovycsv:groovycsv:1.1')
import static com.xlson.groovycsv.CsvParser.parseCsv
import java.io.File

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                Minimap2 Alignment Genome                                                                //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//Data alignment is performed serial with all new files an iteration of file data_aignment_prep is extracting

process minimap_alignment{
	memory "8 GB"
    
    input:
	val string
	val string_array
	
    output:
	val bam_string, emit: bam_string
	val string_array, emit: string_array
    
   
	
	script:
	bam_string=""
    println "Minimap2 Alignment started"
    for(i in 0..string_array.size()-1){
       if (string_array.size() > 0){
           element = string_array.get(i)
           if (params.barcoded == 1){
           basis = element.getParent().getName()
           }
           else {
               basis = element.getParent().getParent().getParent().getName()
           }
           filename = element.getName()
           converted_filename =  filename.minus(".${params.suffix}")
           outname=params.output_dir + basis + "/bam_files/" + converted_filename + ".bam" + " "
           bam_string=bam_string + outname
       }
    }
    data_string = ""
    for (i in 0..params.sample_names.size()-1){
	    data_string = data_string + params.output_dir + params.sample_names.get(i) + "/bam_files/" + " "
    }
    if (string != "")
	    """
        start=\$(date +%s)
    	for i in ${string}
    	do 
            if [ ${params.barcoded} -eq 1 ]
            then
                basis=\$(basename \$(dirname "\${i}"))
                filename=\$(basename "\${i}")
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.output_dir}\${basis}/bam_files/\${converted_filename}
            else
                basis=\$(basename \$(dirname \$(dirname \$(dirname "\${i}"))))
                filename=\$(basename "\${i}")
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.output_dir}\${basis}/bam_files/\${converted_filename}
            fi
            if [ ${params.drs} -eq 1 ]
            then
            minimap2 --MD -ax splice -uf -k14 -t ${params.threads} ${params.output_dir}MT-human_ont.mmi \$i | samtools view -hbS -F 3844 | samtools sort > \${outname} || echo "\$i " >> ${params.output_dir}error_logs/minimap2_genome_failed.log
            else
            minimap2 --MD -ax splice -t ${params.threads} ${params.output_dir}MT-human_ont.mmi \$i | samtools view -hbS -F 3844 | samtools sort > \${outname} || echo "\$i " >> ${params.output_dir}error_logs/minimap2_genome_failed.log 
            fi
        done
        end=\$(date +%s)
        time=\$(echo "\$((\$end-\$start))")
        echo "minimap_genome,${iteration.value},\$time" >> ${params.output_dir}processing_time_table.csv
        function samtoolsParallel  {
        if [ ! \$(ls \${1}| wc -l) -eq 0 ]
        then
            sample=\$1
            run_dir=\$2
            par_basis=\$(basename \$(dirname \${sample}))
            if [ ! -f \${run_dir}bam_genome_merged/\${par_basis}.bam ]
            then
                samtools merge \${run_dir}\${par_basis}/all_gene.bam \${sample}*.bam -f --threads ${params.threads} -c -p || echo "\${sample}" >> \${run_dir}error_logs/merge_genome_all_error.log
                cp \${run_dir}\${par_basis}/all_gene.bam \${run_dir}bam_genome_merged/\${par_basis}.bam 
                samtools index \${run_dir}bam_genome_merged/\${par_basis}.bam || echo "Indexing failed" >> \${run_dir}/error_logs/indexing_genome_all_error.log
            else
            bam_files_to_merge=\${run_dir}\${par_basis}/all_gene.bam" "
            for i in ${bam_string}
            do
                if [[ "\$(basename \$(dirname \$(dirname "\${i}")))" == "\${par_basis}" ]]
                then
                    bam_files_to_merge=\${bam_files_to_merge}"\${i}"" "
                fi
            done
            samtools merge \${run_dir}bam_genome_merged/\${par_basis}.bam \${bam_files_to_merge} -f --threads ${params.threads} -c -p || echo "\${bam_files_to_merge}" >> \${run_dir}error_logs/merge_genome_few_error.log
            cp \${run_dir}bam_genome_merged/\${par_basis}.bam \${run_dir}\${par_basis}/all_gene.bam
            samtools index \${run_dir}bam_genome_merged/\${par_basis}.bam || echo "Indexing failed" >> \${run_dir}/error_logs/indexing_genome_few_error.log
            fi 
        fi
        }
        export -f samtoolsParallel 
        start=\$(date +%s)
        parallel -v -u --env samtoolsParallel --no-notice -j ${params.threads} samtoolsParallel ::: ${data_string} ::: ${params.output_dir}
        end=\$(date +%s)
        time=\$(echo "\$((\$end-\$start))")
        echo "samtools_genome,${iteration.value},\$time" >> ${params.output_dir}processing_time_table.csv
        echo -e "Sample\tnum_reads\tnum_mapped_reads" > ${params.output_dir}mapping_stats.txt
        for i in \$(ls ${params.output_dir}bam_genome_merged/*.bam)
        do
            echo \$i
            numReads=\$(samtools view -c \$i)
            numMappedReads=\$(samtools view -c -F 260 \$i)
            sample=\$(basename \$i) 
            sample=\${sample/.bam/}
            echo "\$sample\t\$numReads\t\$numMappedReads" >> ${params.output_dir}mapping_stats.txt 
        done
        if [ \$(basename \$(dirname \$(dirname \$(dirname${params.script_dir}))))=="NanopoReaTA" ]
        then
            if [ \$(basename \$(dirname \$(dirname${params.script_dir})))=="app" ]
            then
            rm ${params.script_dir}work/* -r || echo "Erase data in nanopore script directory failed"
            rm ${params.script_dir}.nextflow* -r || echo "Erase nextflow log in script directory failed"
            rm ${params.script_dir}../../work/* -r || echo "Erase data in nanopore work directory failed"
            rm ${params.script_dir}../../.nextflow* -r || echo "Erase nextflow log failed"
            fi
        fi
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
//                                                            Minimap2 Alignment Transcriptome                                                             //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process minimap_transcript_alignment{
    memory "8 GB"
	input:
	val string 
	val string_array 	
	output:
	val bam_string, emit: bam_string 
	val string_array, emit: string_array 
	
	script:
    println "Minimap2 Alignment started"
    bam_string=""
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
           outname=params.output_dir + basis + "/bam_files_transcripts/" + converted_filename + ".bam" + " "
           bam_string=bam_string + outname
       }
    }

	if (string != "")
        """
        start=\$(date +%s)
    	for i in ${string}
    	do 
            if [ ${params.barcoded} == 1 ]
            then
                basis=\$(basename \$(dirname \${i}))
                filename=\$(basename \${i})
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.output_dir}\${basis}/bam_files_transcripts/\${converted_filename}
            else
                basis=\$(basename \$(dirname \$(dirname \$(dirname \${i}))))
                filename=\$(basename \${i})
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.output_dir}\${basis}/bam_files_transcripts/\${converted_filename}
            fi
            if [ ${params.drs} -eq 1 ]
            then
            minimap2 --MD -ax map-ont -uf -k14 -t ${params.threads} ${params.output_dir}MT-human_transcript_ont.mmi ${params.output_dir} \$i | samtools view -hbS -F 3844 | samtools sort > \${outname} || echo "\$i " >> ${params.output_dir}error_logs/minimap2_transcript_failed.log
            else
            minimap2 --MD -ax map-ont -t ${params.threads} ${params.output_dir}MT-human_transcript_ont.mmi ${params.output_dir} \$i | samtools view -hbS -F 3844 | samtools sort > \${outname} || echo "\$i " >> ${params.output_dir}error_logs/minimap2_transcript_failed.log
            fi
        done
        end=\$(date +%s)
        time=\$(echo "\$((\$end-\$start))")
        echo "minimap_transcriptome,${iteration.value},\$time" >> ${params.output_dir}processing_time_table.csv
        """
    else
        """
        echo ""
        """
}