version 1.0


# Copies fields from the truvari collapsed inter-sample VCF to the re-genotyped
# inter-sample VCF.
# 
# Remark: the program assumes that CHROM,POS,REF,ALT is a unique identifier of
# each record.
#
workflow SV_Integration_TransferIntersampleAnnotations {
    input {
        File intersample_truvari_bcf
        File intersample_truvari_csi
        File intersample_kanpig_bcf
        File intersample_kanpig_csi
        String columns_to_transfer = 'FORMAT/SUPP_PBSV,FORMAT/SUPP_SNIFFLES,FORMAT/SUPP_PAV,FORMAT/SCORE,FORMAT/CALIBRATION_SENSITIVITY'

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            intersample_truvari_bcf = intersample_truvari_bcf,
            intersample_truvari_csi = intersample_truvari_csi,
            intersample_kanpig_bcf = intersample_kanpig_bcf,
            intersample_kanpig_csi = intersample_kanpig_csi,
            columns_to_transfer = columns_to_transfer,
            docker_image = docker_image
    }
    
    output {
        File output_bcf = Impl.output_bcf
        File output_csi = Impl.output_csi
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
        File intersample_kanpig_bcf
        File intersample_kanpig_csi
        String columns_to_transfer

        Int n_cpu = 2
        Int ram_size_gb = 4
        String docker_image
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 4*ceil( size(intersample_truvari_bcf, "GB") + size(intersample_kanpig_bcf, "GB") )
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations ~{intersample_truvari_bcf} --columns ~{columns_to_transfer} --output-type b --write-index ~{intersample_kanpig_bcf} --output out.bcf
    >>>
    
    output {
        File output_bcf = "out.bcf"
        File output_csi = "out.bcf.csi"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
