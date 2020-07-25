# File Name: Mutect2-TargetCapture-M-PrepareBams.wdl
# Created By: ZWeber
# Created On: 2020-07-22
# Purpose: variant calling on Hybrid Capture Sequencing  Using multisample feature of M2
#   This file is a sub-workflow module meant to take an array of BAMS -> subset them,
#   reheader them, and index them for use in variant calling and other interval based
#   genomics tasks.
# LICENSE: software not developed by Z.Weber retains original licenses.
#   while software developed by Z.Weber can be used under the conditions
#   outlined by the GNU-GPLv3.0

# Workflow Definition

workflow Prepare_Bams {
    # Global Inputs
    Array[File] Bams # array of bam files to process, passed from main pipeline
    Array[File] SampleIDs # array of sample IDs corresponding to BAMS for consistent naming
    File Intervals # intervals over which to subset the bam files
    String WorkflowDocker # name of docker containing gatk and samtools
    Int WorkflowMemory # integer number of GB memory to request on compute node
    Int WorkflowDisk # integer number of GB disk space to requiest on compute node
    Int WorkflowCPU

    # scatter each BAM file onto its own process
    scatter (BamInfo in zip(Bams, SampleIDs)) {
        call Prepare_Bam {
            input:
                bam_file=BamInfo.left,
                bam_name=BamInfo.right,
                sub_targets=Intervals,
                docker_name=WorkflowDocker,
                task_memory=WorkflowMemory,
                task_disk=WorkflowDisk,
                task_cpu=WorkflowCPU
        }
    }

    # output prepped bam, prepped bam index, and the file manifest
    output {
        Array[File] PreppedBam = Prepare_Bam.PreppedBam
        Array[File] PreppedBamIndex = Prepare_Bam.PreppedBamIndex
        Array[File] FileManifest = Prepare_Bam.FileManifest
    }
}


# Task Definitions

## task name: Prepare_Bams
## purpose: miscellaneous tasks related to preparing BAMs for parallel
##      variant calling. Subset bam by target intervals, reheader new
##      bam with apprpriate sample name for consistency, and finally
##      generates a bam index file
task Prepare_Bam {
    #Local Inputs
    File bam_file # bam file to prepare (<b1.bam>)
    String bam_name # sample name of above bam file (<b1name>)
    File sub_targets # # target intervals on which to subset (<sub-name>.interval_list)
    String docker_name # name of docker image containing gatk and samtools
    Int task_memory # GB of task memory required
    Int task_disk # GB of task disk required
    Int task_cpu # integer number of compute cores to reserve for this task
    Int add_cpu = task_cpu - 1 # need the number of cores minus 1 for samtools

    # script execution
    command {
        # convert target interval_list format into BED format for samtools
        gatk --java-options "-Xmx${task_memory}G" IntervalListToBed \
            -I ${sub_targets} \
            -O "sub_targets.bed"

        # use samtools to reheader, subset and index the passed bam
        samtools view -@ ${add_cpu} -H ${bam_file} \
            | sed -e 's/\(SM:\)\+\(\S\)\+/SM:'${bam_name}'/g' \
            | samtools reheader - ${bam_file} \
            | samtools view -@ ${add_cpu} -L "sub_targets.bed" -b - \
            > ${bam_name}'-fixed-subset.bam'
        samtools index -@ ${add_cpu} -b ${bam_name}'-fixed-subset.bam'

        # for debugging purposes, generate a file manifest
        ls -lh > "MANIFEST.txt"
    }

    # specify runtime params
    runtime {
        docker: docker_name
        memory: task_memory + "G"
        disks: "local-disk " + task_disk + " HDD"
        cpu: task_cpu
    }

    # output the prepared bam and its index, as well as the file manifest
    output {
        File PreppedBam = bam_name + '-fixed-subset.bam'
        File PreppedBamIndex = bam_name + '-fixed-subset.bam.bai'
        File FileManifest = "MANIFEST.txt"
    }
}
