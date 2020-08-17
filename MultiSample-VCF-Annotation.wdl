# File Name: MultiSample-VEP-Annotation.wdl
# Created By: ZWeber
# Created On: 2020-05-31
# Purpose: Annotate Position-based variants from a Multisample VCF
# LICENSE: software not developed by Z.Weber retains original licenses.
#   while software developed by Z.Weber can be used under the conditions
#   outlined by the GNU-GPLv3.0

workflow multisample_vcf_annotate {
    # Global Inputs
    File MultiSample_Vcf # multisample vcf to annotate
    File MultiSample_Vcf_Index # index file of vcf to annotate
    File VEP_Anno_DB # annotation resource for VEP (should be Ensembl resource)
    String Sample_Set_ID # Identifier to serve as prefix for generated files
    String VEP_Docker # name of docker image for this workflow
    Int BootDisk_Size # amount of GB to give as additional boot size
    Int Memory_Size #  amount of GB to give memory
    Int Disk_Size # amount of Hard Disk space to allocate

    # call annotation task
    call annotate_variants {
        input:
            vcf=MultiSample_Vcf,
            vcf_index=MultiSample_Vcf_Index,
            vep_datasource=VEP_Anno_DB,
            sample_set_id=Sample_Set_ID,
            docker_name=VEP_Docker,
            boot_disk=BootDisk_Size,
            task_memory=Memory_Size,
            task_disk=Disk_Size
    }

    # Global Outputs (annotated vcf and annotated table)
    output {
        File VEP_Annotated_VCF = annotate_variants.annotated_vcf
        File VEP_Annotated_VCF_Idx = annotate_variants.annotated_vcf_idx
        File VEP_Annotated_Table = annotate_variants.annotated_table
        File VEP_Annotated_Table_DataDict = annotate_variants.annotated_table_datadict
    }
}

# task annotate_variants takes a multisample VCF and
# annoates using VEP before returning a number of files
task annotate_variants {
    # local inputs
    File vcf # multisample vcf to annoate
    File vcf_index # multisample vcf index
    File vep_datasource # ensembl database (.tar.gz)
    String sample_set_id # file name prefix for sample size
    String docker_name # name of docker image
    Int boot_disk # size of extra boot disk for the process
    Int task_memory # GB of memory required for the task
    Int task_disk # GB of disk required for the task

    # run annotation task
    command {
        # unzip database here in place
        tar -xzf ${vep_datasource}

        # gunzip the vcf if required
        python3 /home/scripts/prepare_multivcf.py \
            -i ${vcf} \
            -o ${sample_set_id}".vcf"

        # run VEP to generate annotated VCF
        vep --offline --dir_cache `pwd` \
            -i ${sample_set_id}".vcf" --format "vcf" \
            -o ${sample_set_id}".anno.vcf" --vcf \
            --canonical --protein --check_existing \
            --symbol --regulatory --pick --biotype \
            --polyphen b --sift b

        # run VEP to generate tabular file
        vep --offline --dir_cache `pwd` \
            -i ${sample_set_id}".vcf" --format "vcf" \
            -o ${sample_set_id}".anno.txt" --tab \
            --canonical --protein --check_existing \
            --symbol --regulatory --pick --biotype \
            --polyphen b --sift b

        # post-process tabular VEP results and VCF annotation
        bcftools query -f "%FILTER\t%INFO/TLOD\t%INFO/NLOD\t%INFO/NALOD[\t%DP\t%AD]\n" \
            ${sample_set_id}".anno.vcf" > "variant-call-info.txt"
        bcftools query -l ${sample_set_id}".anno.vcf" > "sample-names.txt"

        python3 /home/scripts/build-variant-annotation-table.py \
            --tab ${sample_set_id}".anno.txt" \
            --vcf-info "variant-call-info.txt" \
            --sample-names "sample-names.txt" \
            --output-filename-prefix ${sample_set_id}

        # tabix index and gzip anootated vcf before returning
        bgzip ${sample_set_id}".anno.vcf"
        tabix ${sample_set_id}".anno.vcf.gz"
    }

    # runtime (small task with large boot disk)
    runtime {
        docker: docker_name
        memory:  task_memory + "G"
        disks: "local-disk " + task_disk + " HDD"
        bootDiskSizeGb: boot_disk
    }

    # output annotated files
    output {
        File annotated_vcf = sample_set_id + ".anno.vcf.gz"
        File annotated_vcf_idx = sample_set_id + ".anno.vcf.gz.tbi"
        File annotated_table = sample_set_id + "-annotation-table.txt"
        File annotated_table_datadict = sample_set_id + "-datadict.txt"
    }
}
