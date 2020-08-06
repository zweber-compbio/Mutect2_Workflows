# File Name: Mutect2-TargetCapture-Pipeline-VariantCalling.wdl
# Created By: Zachary Weber
# Created On: 2020-07-22
# Last Updated: 2020-08-04
# Purpose: Short variant calling on Hybrid Capture Sequencing data (WES/TPS)
# implementing M2 Multisample Mode with OB filter and force calling options
# LICENSE: software not developed by Z.Weber retains original licenses.
#   while software developed by Z.Weber can be used under the conditions
#   outlined by the GNU-GPLv3.0

# import external modules for sub-workflows
import "https://github.com/zweber-compbio/Mutect2_Workflows/blob/master/Mutect2-TargetCapture-M-PrepareBams.wdl" as Prepare_Bams


# Workflow Definition
workflow Multisample_Variant_Calling {
    # ---------------------
    # Global Inputs
    ## genome reference related
    File ReferenceFasta # reference sequence fasta file
    File ReferenceFastaIndex # reference sequence fasta index from samtools
    File ReferenceFastaDict # reference sequence dictionary from samtools

    ## target capture and parallelization
    File TargetIntervals # .interval_list file containing capture regions
    Int TargetScatterCount # number of parallel partitions of variant caller
    Int TargetIntervalsPadding # amount of extra bp to pad interval shoulders
    String? ExcludeRegions # regions of the genome to ignore becuase they are problematic

    ## bam files and sample Names
    Array[File] TumorBams # array of tumor bams for variant caller
    Array[File] NormalBams # array of matched normals for variant caller
    Array[String] TumorBamSampleNames # array of bam sample names in same order as bams
    Array[String] NormalBamSampleNames # array of bam sample names in same ordre as bams
    String SampleSetName # name of sample set to be used in output file naming

    ## vcf data for PON, germline population frequency ref, force calling alleles
    File Mutect2PON # panel of normals vcf generated by mutect2 on germline
    File Mutect2PONIndex # index file for panel of normals generated by mutect2
    File GnomadVcf # allele frequency only VCF for germline filtering in pop
    File GnomadVcfIndex # allele frequency only VCF index
    File? ForceCallSitesVcf # vcf for known variants to assess (should be present im force-call mode)
    File? ForceCallSitesVcfIndex # vcf index (should be present im force-call mode)

    ## modes of operation for pipeline at large
    Boolean UseForceCalling = false # whether or not to force call passed sites
    Boolean UseOBModel = true # apply orientation bias model filter (useful with ffpe)

    ## runtime parameters for cloud compute
    String VariantCallDocker # docker image playing host to variant calling tools
    Int LargeCPU = 8 # number of CPUS on compute intensive tasks
    Int SmallMemory = 8 # GB of memory allocated to small tasks in workflow
    Int LargeMemory = 16 # GB of memory allocated to large tasks in workflow
    Int SmallDisk = 100 # GB of disk space allocated to small tasks
    Int LargeDisk = 100 # GB of disk space allocated to large tasks

    # ---------------------
    ## 1. split target intervals into subsets based on genome region
    ##    and scatter those intervals into separate processes
    ##    for purposes of parallelization

    call Split_Targets {
        input:
            ref_fa=ReferenceFasta,
            ref_fai=ReferenceFastaIndex,
            ref_dict=ReferenceFastaDict,
            targets=TargetIntervals,
            target_padding=TargetIntervalsPadding,
            exclude_regions=ExcludeRegions,
            scatter_count=TargetScatterCount,
            docker_name=VariantCallDocker,
            task_memory=SmallMemory,
            task_disk=SmallDisk
    }

    scatter (SubTargets in Split_Targets.SplitIntervals) {

        # ---------------------
        ## 2A,B,C. subset PON, germline resource, and force calling regions
        ##         by new, smaller targets set provided in previous step

        call Subset_Vcf as Subset_M2PON {
            input:
                sub_targets=SubTargets,
                output_vcf_name="panel-of-normals",
                vcf_file=Mutect2PON,
                vcf_file_index=Mutect2PONIndex,
                docker_name=VariantCallDocker,
                task_memory=SmallMemory,
                task_disk=SmallDisk
        }

        call Subset_Vcf as Subset_GnomadGermlineResource {
            input:
                sub_targets=SubTargets,
                output_vcf_name="gnomad",
                vcf_file=GnomadVcf,
                vcf_file_index=GnomadVcfIndex,
                docker_name=VariantCallDocker,
                task_memory=SmallMemory,
                task_disk=SmallDisk
        }

        if (UseForceCalling) {
            call Subset_Vcf as Subset_ForceCallSites {
                input:
                    sub_targets=SubTargets,
                    output_vcf_name="force-calls",
                    vcf_file=ForceCallSitesVcf,
                    vcf_file_index=ForceCallSitesVcfIndex,
                    docker_name=VariantCallDocker,
                    task_memory=SmallMemory,
                    task_disk=SmallDisk
            }
        }

        # ---------------------
        ## 2D. subset tumor BAM files by new, smaller
        ##     target set provided in step 1. Also, reheader and
        ##     index these smaller BAM files to ensure consistent
        ##     sample naming across parallel processes

        call Prepare_Bams.Prepare_Bams as Prepare_TumorBams {
            input:
                Bams=TumorBams,
                SampleIDs=TumorBamSampleNames,
                Intervals=SubTargets,
                WorkflowDocker=VariantCallDocker,
                WorkflowMemory=SmallMemory,
                WorkflowDisk=SmallDisk,
                WorkflowCPU=LargeCPU
        }

        # ---------------------
        ## 2E. subset normal BAM files by new, smaller
        ##    target set provided in step 1. Also, reheader and
        ##    index these smaller BAM files to ensure consistent
        ##    sample naming across parallel processes

        call Prepare_Bams.Prepare_Bams as Prepare_NormalBams {
            input:
                Bams=NormalBams,
                SampleIDs=NormalBamSampleNames,
                Intervals=SubTargets,
                WorkflowDocker=VariantCallDocker,
                WorkflowMemory=SmallMemory,
                WorkflowDisk=SmallDisk,
                WorkflowCPU=LargeCPU
        }

        # ---------------------
        ## 3. run mutect2 variant caller on target sets provided in step 1.
        ##    this task may require larger than average memory and cpu for
        ##    efficient run times

        call Mutect2_Call_Variants {
            input:
                ref_fa=ReferenceFasta,
                ref_fai=ReferenceFastaIndex,
                ref_dict=ReferenceFastaDict,
                sub_targets=SubTargets,
                tumor_bams=Prepare_TumorBams.PreppedBam,
                tumor_bam_indices=Prepare_TumorBams.PreppedBamIndex,
                normal_bams=Prepare_NormalBams.PreppedBam,
                normal_bam_indices=Prepare_NormalBams.PreppedBamIndex,
                normal_bam_names=NormalBamSampleNames,
                mutect2_pon=Subset_M2PON.SubsetVcf,
                mutect2_pon_index=Subset_M2PON.SubsetVcfIndex,
                gnomad_vcf=Subset_GnomadGermlineResource.SubsetVcf,
                gnomad_vcf_index=Subset_GnomadGermlineResource.SubsetVcfIndex,
                force_calling_vcf=Subset_ForceCallSites.SubsetVcf,
                force_calling_vcf_index=Subset_ForceCallSites.SubsetVcfIndex,
                docker_name=VariantCallDocker,
                task_memory=LargeMemory,
                task_disk=LargeDisk,
                task_cpu=LargeCPU
        }
    }

    # ---------------------
    ## 4. merge the unfiltered mutect calls, and call stat files from
    ##    all of the scattered variant calling tasks

    call Merge_Mutect2_Results {
        input:
            unfiltered_call_vcfs=Mutect2_Call_Variants.Mutect2UnfiltVcf,
            unfiltered_call_vcf_indices=Mutect2_Call_Variants.Mutect2UnfiltVcfIndex,
            call_stat_files=Mutect2_Call_Variants.Mutect2CallStats,
            docker_name=VariantCallDocker,
            task_memory=SmallMemory,
            task_disk=SmallDisk
    }

    # ---------------------
    ## 5. If we have the orientation bias filter turned on, take the
    ##    OB metrics collected during step 3, and learn an orientation
    ##    bias model for downstream filtering tasks

    if (UseOBModel) {
        call Learn_ReadOBModel {
            input:
                orientation_bias_files=Mutect2_Call_Variants.Mutect2F1R2Metrics,
                docker_name=VariantCallDocker,
                task_memory=SmallMemory,
                task_disk=SmallDisk
        }
    }

    # ---------------------
    ## 6. Filter the raw mutect2 calls to label false positives,
    ##    germline population snps, PON artifacts, OB artifacts, etc.

    call Filter_Mutect2_Calls {
        input:
            ref_fa=ReferenceFasta,
            ref_fai=ReferenceFastaIndex,
            ref_dict=ReferenceFastaDict,
            sample_set_name=SampleSetName,
            unfiltered_vcf=Merge_Mutect2_Results.MergedUnfiltVcf,
            unfiltered_vcf_index=Merge_Mutect2_Results.MergedUnfiltVcfIndex,
            ob_prior_table=Learn_ReadOBModel.OBPriors,
            call_stats_file=Merge_Mutect2_Results.MergedMutect2CallStats,
            docker_name=VariantCallDocker,
            task_memory=SmallMemory,
            task_disk=SmallDisk
    }

    # ---------------------
    ## 7. Collect Force Called Alleles
    ##    If we used force-calling mode, create a vcf file that contains
    ##    the results of calls at those loci

    if (UseForceCalling) {
        call Get_ForceCalled_Alleles {
            input:
                force_calling_vcf=ForceCallSitesVcf,
                force_calling_vcf_index=ForceCallSitesVcfIndex,
                filter_annotated_vcf=Filter_Mutect2_Calls.FilteredMutect2Calls,
                filter_annotated_vcf_index=Filter_Mutect2_Calls.FilteredMutect2CallsIndex,
                sample_set_name=SampleSetName,
                docker_name=VariantCallDocker,
                task_memory=SmallMemory,
                task_disk=SmallDisk
        }
    }

    # ---------------------
    ## 8. Output Final Files
    ## Output Globally:
    ##      * Unfiltered Vcf file plus index
    ##      * Filtered Vcf file plus index
    ##      * Passing Only Vcf file plus index
    ##      * Force Calling Vcf file plus index

    output {
        File Mutect2_Unfiltered_Calls = Merge_Mutect2_Results.MergedUnfiltVcf
        File Mutect2_Unfiltered_Calls_Index = Merge_Mutect2_Results.MergedUnfiltVcfIndex
        File Mutect2_Filtered_Calls = Filter_Mutect2_Calls.FilteredMutect2Calls
        File Mutect2_Filtered_Calls_Index = Filter_Mutect2_Calls.FilteredMutect2Calls
        File Mutect2_Passing_Calls = Filter_Mutect2_Calls.PassingMutect2Calls
        File Mutect2_Passing_Calls_Index = Filter_Mutect2_Calls.PassingMutect2CallsIndex
        File? Mutect2_ForceCalled_Alleles = Get_ForceCalled_Alleles.FilteredForceCalledAlleles
        File? Mutect2_ForceCalled_Alleles_Index = Get_ForceCalled_Alleles.FilteredForceCalledAllelesIndex
    }
}


#  Task Definitions

## task name: Split_Targets
## purpose: divides task into smaller genome regions for
##     efficient and cost effective parallel processing
task Split_Targets {
    # Local Inputs
    File ref_fa # reference sequence fasta file (<ref-name>.fasta)
    File ref_fai # reference sequence fasta index from samtools (<ref-name>.fasta.fai)
    File ref_dict # reference sequence dictionary from samtools (<ref-name>.dict)
    File targets # file containing sequencing regions of interest (<regions>.interval_list)
    String? exclude_regions # intervals to exclude from the master list before splitting
    Int target_padding # integer number of bp to pad targets on each side
    Int scatter_count # number of sub regions (and subsequent parallel processes) to divide regions
    String docker_name # name of docker image contianing gatk toolchain
    Int task_memory # number of GB of memory for this task
    Int task_disk # number of GB of disk space for this task

    # script execution
    command {
        mkdir subintervals/
        gatk --java-options "-Xmx${task_memory}G" "SplitIntervals" \
            --output subintervals \
            --reference ${ref_fa} \
            --intervals ${targets} \
            --interval-padding ${target_padding} \
            --scatter-count ${scatter_count} \
            ${default="" "-XL " + exclude_regions}
    }

    # specify runtime params
    runtime {
        docker: docker_name
        memory: task_memory + "G"
        disks: "local-disk " + task_disk + " HDD"
    }

    # output array of split target intervals
    output {
        Array[File] SplitIntervals = glob("subintervals/*.interval_list")
    }
}


## task name: Subset_Vcf
## purpose: subset a vcf of interest by the provided calling region
##      to reduce disk footprint and time spent localizing files to nodes
task Subset_Vcf{
    # Local Inputs
    File sub_targets # target intervals on which to subset (<sub-name>.interval_list)
    File vcf_file # vcf file to subset based on genomic intervals
    File vcf_file_index # # index file to the vcf of interest
    String output_vcf_name # string representing file prefix for ouput vcf
    String docker_name # name of docker image contianing gatk toolchain
    Int task_memory # GB of task memory required
    Int task_disk # GB of task disk required

    # script execution
    command {
        # run gatk SelectVariants on vcf of interest
        gatk --java-options "-Xmx${task_memory}G" "SelectVariants" \
            -V ${vcf_file} \
            -L ${sub_targets} \
            -O ${output_vcf_name}"-subset.vcf.gz"
    }

    # specify runtime params
    runtime {
        docker: docker_name
        memory: task_memory + "G"
        disks: "local-disk " + task_disk + " HDD"
    }

    # output the subset vcf file and its corresponding index
    output {
        File SubsetVcf = output_vcf_name + "-subset.vcf.gz"
        File SubsetVcfIndex = output_vcf_name + "-subset.vcf.gz.tbi"
    }
}


## task name: Mutect2_Call_Variants
## purpose: Calls somatic variants using mutect2 on the regions
## provided in the scattered intervals. these calls are unfiltered
task Mutect2_Call_Variants {
    # Local Inputs
    File ref_fa # reference sequence fasta file (<ref-name>.fasta)
    File ref_fai # reference sequence fasta index from samtools (<ref-name>.fasta.fai)
    File ref_dict # reference sequence dictionary from samtools (<ref-name>.dict)
    File sub_targets # file containing sequencing regions over which to call (<regions>.interval_list)
    Array[File] tumor_bams # array of prepped bams [<b1-prepped>.bam, <b2-prepped>.bam, ...]
    Array[File] tumor_bam_indices # array of index files to prepped bams [<b1-prepped>.bam.bai, ...]
    Array[File] normal_bams # array of prepped bams [<b1-prepped>.bam, <b2-prepped>.bam, ...]
    Array[File] normal_bam_indices # array of index files to prepped bams [<b1-prepped>.bam.bai, ...]
    Array[String] normal_bam_names # array of strings matching normal bam sample IDS
    File mutect2_pon # panel of normals vcf generated on germline samples
    File mutect2_pon_index # index file for mutect2_pon
    File gnomad_vcf # gnomad allele frequency resource
    File gnomad_vcf_index # gnomad allele frequency resource index file
    File? force_calling_vcf # vcf for sites to force call if defined in this workflow
    File? force_calling_vcf_index # index to the above vcf if defined in this workflow

    String docker_name # name of docker image containing gatk
    Int task_memory # GB of task memory required
    Int task_disk # GB of task disk required
    Int task_cpu # integer number of compute cores to reserve for this task

    # script execution
    command {
        # run mutect2 on adjusted BAM files
        gatk --java-options "-Xmx${task_memory}G" Mutect2 \
            -O "somatic.vcf.gz" \
            --reference ${ref_fa} \
            -I ${sep=" -I " tumor_bams} \
            -I ${sep=" -I " normal_bams} \
            -normal ${sep=" -normal " normal_bam_names} \
            -L ${sub_targets} \
            -pon ${mutect2_pon} \
            --germline-resource ${gnomad_vcf} \
            --genotype-pon-sites \
            --f1r2-tar-gz "f1r2.tar.gz" \
            --genotype-germline-sites ${default="" "--alleles " + force_calling_vcf}
    }

    # specify runtime params
    runtime {
        docker: docker_name
        memory: task_memory + "G"
        disks: "local-disk " + task_disk + " HDD"
        cpu: task_cpu
    }

    # output somatic variants plus orientation bias model information
    output {
        File Mutect2UnfiltVcf="somatic.vcf.gz"
        File Mutect2UnfiltVcfIndex="somatic.vcf.gz.tbi"
        File Mutect2CallStats="somatic.vcf.gz.stats"
        File Mutect2F1R2Metrics="f1r2.tar.gz"
    }
}


## task name: Merge_Mutect2_Results
## purpose: take unfiltered calls and call stats files from scattered
##      variant calling processes and merge into a single indexed vcf
##      as well as a single unified call stats file.
task Merge_Mutect2_Results {
    # Local Inputs
    Array[File] unfiltered_call_vcfs # unfiltered vcf files produced by scattered variant calling
    Array[File] unfiltered_call_vcf_indices # index files corresponding to unfiltered vcfs
    Array[File] call_stat_files # call stats output produced by mutect2
    String docker_name # name of docker image containing gatk
    Int task_memory # GB of task memory required
    Int task_disk # GB of task disk required

    # script execution
    command {
        # merge calls
        gatk --java-options "-Xmx${task_memory}G" MergeVcfs \
            -I ${sep=' -I ' unfiltered_call_vcfs} \
            -O "merged.somatic.vcf.gz"

        # merge call stats
        gatk --java-options "-Xmx${task_memory}G" MergeMutectStats \
            -stats ${sep=" -stats " call_stat_files} \
            -O "merged.mutect2_calls.stats"
    }
    # specify runtime params
    runtime {
        docker: docker_name
        memory: task_memory + "G"
        disks: "local-disk " + task_disk + " HDD"
    }

    # ouput merged unfiltered calls and call stat files
    output {
        # output merged unfiltered calls as well as call stats
        File MergedUnfiltVcf = "merged.somatic.vcf.gz"
        File MergedUnfiltVcfIndex = "merged.somatic.vcf.gz.tbi"
        File MergedMutect2CallStats = "merged.mutect2_calls.stats"
    }
}


## task name: Learn_ReadOBModel
## purpose: takes context specific information from call sites and
##      generates a model used for filtering chemical artifacts
##      introduced in  library preparation, such as those caused
##      by FFPE storage
task Learn_ReadOBModel {
    # local linputs
    Array[File] orientation_bias_files # array of f1r1 files generated during Mutect2 calling
    String docker_name # name of docker image containing gatk
    Int task_memory # GB of task memory required
    Int task_disk # GB of task disk required

    # script execution
    command {
        # learn ob model
        gatk --java-options "-Xmx${task_memory}G" LearnReadOrientationModel \
                -I ${sep=" -I " orientation_bias_files} \
                -O "artifact-priors.tar.gz"
    }

    # run using small disk specification
    runtime {
        docker: docker_name
        memory: task_memory + "G"
        disks: "local-disk " + task_disk + " HDD"
    }

    # output the artifact priors table
    output {
        File OBPriors = "artifact-priors.tar.gz"
    }
}


## task name: Filter_Mutect2_Calls
## purpose: Filter false positives, germline variants, and other
##      seuqencing artifacts from mutect2 outputs, as well as create
##      an output vcf for only those variants which pass all standard
##      filters
task Filter_Mutect2_Calls {
    # Local Inputs
    File ref_fa # reference sequence fasta file (<ref-name>.fasta)
    File ref_fai # reference sequence fasta index from samtools (<ref-name>.fasta.fai)
    File ref_dict # reference sequence dictionary from samtools (<ref-name>.dict)
    String sample_set_name # string representing the name of the sample set, for output file naming
    File unfiltered_vcf # unfiltered vcf file from mutect2, merged across parallel callsets
    File unfiltered_vcf_index # index file to the merged mutect2 vcf output
    File? ob_prior_table # priors table from the learn OB model task
    File call_stats_file # merged call stats file from mutect2
    String docker_name # name of docker image containing gatk
    Int task_memory # GB of task memory required
    Int task_disk # GB of task disk required

    # script execution
    command {
        # filter variants called by mutect2
        gatk --java-options "-Xmx${task_memory}G" FilterMutectCalls \
            --reference ${ref_fa} \
            -V ${unfiltered_vcf} \
            -stats ${call_stats_file} \
            --filtering-stats ${sample_set_name}"-somatic-filtering.stats" \
            -O ${sample_set_name}"-somatic-filt.vcf.gz" \
            ${default="" "-ob-priors " + ob_prior_table}

        # create file of passing variants only
        gatk --java-options "-Xmx${task_memory}G" SelectVariants \
            -V ${sample_set_name}"-somatic-filt.vcf.gz" \
            -O ${sample_set_name}"-somatic-passing.vcf.gz" \
            --exclude-filtered
    }

    # specify runtime params
    runtime {
        docker: docker_name
        memory: task_memory + "G"
        disks: "local-disk " + task_disk + " HDD"
    }

    # output all variants with filtering determinations
    # plus a passing-only vcf
    output {
        File FilteredMutect2Calls = sample_set_name + "-somatic-filt.vcf.gz"
        File FilteredMutect2CallsIndex = sample_set_name + "-somatic-filt.vcf.gz.tbi"
        File PassingMutect2Calls = sample_set_name + "-somatic-passing.vcf.gz"
        File PassingMutect2CallsIndex = sample_set_name + "-somatic-passing.vcf.gz.tbi"
    }
}


## task name: Get_ForceCalled_Alleles
## purpose: extract force calls from the filtered variant set
##      and return them in a separete vcf, regardless of their
##      filter status.
task Get_ForceCalled_Alleles {
    # Local Inputs
    File force_calling_vcf # allleles of interest, force called during mutect2 task
    File force_calling_vcf_index # index file to the above vcf
    File filter_annotated_vcf # merged vcf with filtering annotations
    File filter_annotated_vcf_index # index file to merged, filtered variants
    String sample_set_name # string representing the name of the sample set, for output file naming
    String docker_name # name of docker image containing gatk
    Int task_memory # GB of task memory required
    Int task_disk # GB of task disk required

    # script execution
    command {
        gatk --java-options "-Xmx${task_memory}G" SelectVariants \
            -V ${filter_annotated_vcf} \
            -L ${force_calling_vcf} \
            -O ${sample_set_name}"-force-calls.vcf.gz"
    }

    # specify runtime params
    runtime {
        docker: docker_name
        memory: task_memory + "G"
        disks: "local-disk " + task_disk + " HDD"
    }

    # output the force called alleles and index file
    output {
        File FilteredForceCalledAlleles = sample_set_name + "-force-calls.vcf.gz"
        File FilteredForceCalledAllelesIndex = sample_set_name + "-force-calls.vcf.gz.tbi"
    }
}
