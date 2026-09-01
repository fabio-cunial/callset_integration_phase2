version 1.0


# Subsets a multi-sample VCF to a given set of samples, removing also all
# records that are REF or missing after the subset operation.
#
workflow SubsetVcfToGivenSamples {
    input {
        String prefix
        File input_bcf
        File samples_file

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            prefix = prefix,
            input_bcf = input_bcf,
            samples_file = samples_file,

            docker_image = docker_image
    }
    
    output {
        File output_step1_bcf = Impl.output_step1_bcf
        File output_step1_csi = Impl.output_step1_csi
        File output_step2_bcf = Impl.output_step2_bcf
        File output_step2_csi = Impl.output_step2_csi
    }
}


# Performance on the main VCF (13GB) on a 4-cores HDD:
#
# TASK                                      CPU%       RAM     TIME
# bcftools view --samples-file              200%       500M    1h40m
# bcftools +fill-tags | bcftools view                          1h
# bcftools index                            130%       30M     15m
#
task Impl {
    input {
        String prefix
        File input_bcf
        File samples_file

        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 4
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_bcf, "GB"))
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"

        # Subsetting
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file ~{samples_file} --output-type b1 --write-index=csi ~{input_bcf} --output subset_to_samples_step1.bcf
        date 1>&2
        bcftools +fill-tags --threads ${N_THREADS} --output-type u subset_to_samples_step1.bcf -- --tags AC | bcftools view --threads ${N_THREADS} --min-ac 1 --output-type b --output subset_to_samples_step2.bcf
        date 1>&2
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} subset_to_samples_step2.bcf

        # Outputting
        mv subset_to_samples_step1.bcf ~{prefix}_subset_to_samples_step1.bcf
        mv subset_to_samples_step1.bcf.csi ~{prefix}_subset_to_samples_step1.bcf.csi
        mv subset_to_samples_step2.bcf ~{prefix}_subset_to_samples_step2.bcf
        mv subset_to_samples_step2.bcf.csi ~{prefix}_subset_to_samples_step2.bcf.csi
    >>>
    
    output {
        File output_step1_bcf = prefix + "_subset_to_samples_step1.bcf"
        File output_step1_csi = prefix + "_subset_to_samples_step1.bcf.csi"
        File output_step2_bcf = prefix + "_subset_to_samples_step2.bcf"
        File output_step2_csi = prefix + "_subset_to_samples_step2.bcf.csi"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
