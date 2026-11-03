version 1.0


# Prints histograms of `SUPP_*` and `SCORE` values of a VCF record over all of
# its samples. See `GetDiscoveryAnnotations.java` for details.
#
workflow SV_Integration_GetDiscoveryAnnotations {
    input {
        File intersample_truvari_bcf
        File intersample_truvari_csi

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            intersample_truvari_bcf = intersample_truvari_bcf,
            intersample_truvari_csi = intersample_truvari_csi,
            docker_image = docker_image
    }
    
    output {
        File output_tsv = Impl.output_tsv
    }
}


# Performance with 4 cores and 8 GB of RAM:
#
# TASK                      % CPU       RAM     TIME
# 
#
task Impl {
    input {
        File intersample_truvari_bcf
        File intersample_truvari_csi
        Int n_score_quanta = 6
        Float score_min = 0.4
        Float score_max = 1.0

        Int n_cpu = 2
        Int ram_size_gb = 4
        String docker_image
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 2*ceil( size(intersample_truvari_bcf,"GB") )
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        date 1>&2
        bcftools query --format '%ID[\t%SUPP_PBSV,%SUPP_SNIFFLES,%SUPP_PAV,%SCORE]\n' ~{intersample_truvari_bcf} | java -cp ~{docker_dir} GetDiscoveryAnnotations - ~{n_score_quanta} ~{score_min} ~{score_max} out.tsv
        date 1>&2
    >>>
    
    output {
        File output_tsv = "out.tsv"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
