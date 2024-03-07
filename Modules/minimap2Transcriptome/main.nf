params.fastq_path="/home/vincent/projects/hackathon/HEK_HELA_fastq/HEK_HeLa/20230323_1311_MN32609_FAQ51879_03abc106/fastq_pass/barcode01/FAQ51879_pass_barcode01_03abc106_254e4132_14.fastq.gz"
params.ref_path="/home/vincent//projects/hackathon/data_hackathon/chr20.fa"
params.output_dir="/home/vincent/projects/hackathon"
params.threads=5
params.drs=1

process minimap_transcript_alignment {
    memory "8 GB"
	input:
        val string from minimap_transcript_channel
        val string_array from minimap_transcript_channel2	
	output:
        val bam_string into move_transcripts_channel
        val string_array into move_transcripts_channel2
	
	script:
        println "Minimap2 Alignment started"
        bam_string=""
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
            outname=params.run_dir + basis + "/bam_files_transcripts/" + converted_filename + ".bam" + " "
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
                outname=${params.run_dir}\${basis}/bam_files_transcripts/\${converted_filename}
            else
                basis=\$(basename \$(dirname \$(dirname \$(dirname \${i}))))
                filename=\$(basename \${i})
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.run_dir}\${basis}/bam_files_transcripts/\${converted_filename}
            fi
            if [ ${params.DRS} -eq 1 ]
            then
            minimap2 --MD -ax map-ont -uf -k14 -t ${params.threads} ${params.run_dir}MT-human_transcript_ont.mmi ${params.run_dir} \$i | samtools view -hbS -F 3844 | samtools sort > \${outname} || echo "\$i " >> ${params.run_dir}error_logs/minimap2_transcript_failed.log
            else
            minimap2 --MD -ax map-ont -t ${params.threads} ${params.run_dir}MT-human_transcript_ont.mmi ${params.run_dir} \$i | samtools view -hbS -F 3844 | samtools sort > \${outname} || echo "\$i " >> ${params.run_dir}error_logs/minimap2_transcript_failed.log
            fi
        done
        end=\$(date +%s)
        time=\$(echo "\$((\$end-\$start))")
        echo "minimap_transcriptome,${iteration.value},\$time" >> ${params.run_dir}processing_time_table.csv
        """
    else
        """
        echo ""
        """
}

process move_transcript_files{
    memory "4 GB"
    input: 
    val bam_string from move_transcripts_channel
	val string_array from move_transcripts_channel2

    output:
    val bam_full_string into moving_done_channel
    val string_array into moving_done_channel2

    script:
    bam_full_string = ""
    for(i in 0..string_array.size()-1){
       if (string_array.size() > 0){
           element = string_array.get(i)
           if (params.barcoded == 1){
           basis = element.getParent().getName()
           }
           else {
               basis = element.getParent().getParent().getParent().getName()
           }
           filename= element.getName()

           converted_filename =  filename.minus(".${params.suffix}")
           outname= params.run_dir + basis + "/bam_files_transcripts/" + "full/" + converted_filename
           bam_full_string=bam_full_string + outname + ".bam" + " "
       }
    }
    move_transcripts_minimap2_running.value = 1
    if (bam_string != "")
        """
        for i in ${bam_string}
        do
            basis=\$(dirname \${i})
            filename=\$(basename \${i})
            cp \${i} \$basis/full/\$filename || echo "File does not exist: \$basis/full/\$filename" > ${params.run_dir}error_logs/copy_full_bam_failed.log
        done
        """
    else
        """
        echo empty
        """
}

process moving_done {
     input:
        val bam_full_string from moving_done_channel
        val string_array from moving_done_channel2

     output:
        val bam_full_string into salmon_channel
        val string_array into salmon_channel2

     exec:
        if (bam_full_string != ""){
            move_transcripts_minimap2_running.value = 0
        }
}