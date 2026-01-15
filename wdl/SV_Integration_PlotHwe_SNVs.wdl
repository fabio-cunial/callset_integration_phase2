version 1.0

#
workflow SV_Integration_PlotHwe_SNVs {
    input {
        File intersample_bcf
        File intersample_csi
        
        Float min_af = 0.1
        
        Array[String] ancestry_name
        Array[File] ancestry_samples
        
        File? plothw_r
    }
    parameter_meta {
        ancestry_samples: "A list of sample IDs for each ancestry."
    }
    
    # All
    call SelectBiallelic as all {
        input:
            vcf_gz = intersample_bcf,
            vcf_tbi = intersample_csi
    }
    call Vcf2Counts as all_counts {
        input:
            vcf_gz = all.out_vcf_gz,
            tbi = all.out_tbi
    }
    call Counts2Plot as all_plot {
        input:
            gt_counts = all_counts.gt_counts,
            out_file_name = "all.png",
            plothw_r = plothw_r
    }

    # Frequent
    call FilterByAf as frequent {
        input:
            bcf = all.out_vcf_gz,
            csi = all.out_tbi,
            min_af = min_af,
            smaller_or_larger = 1
    }
    call Vcf2Counts as frequent_counts {
        input:
            vcf_gz = frequent.out_vcf_gz,
            tbi = frequent.out_tbi
    }
    call Counts2Plot as frequent_plot {
        input:
            gt_counts = frequent_counts.gt_counts,
            out_file_name = "frequent.png",
            plothw_r = plothw_r
    }
    
    # All of the above, but by ancestry.
    scatter (i in range(length(ancestry_samples))) {
        # All
        call FilterBySamples as ancestry_all {
            input:
                bcf = all.out_vcf_gz,
                csi = all.out_tbi,
                sample_ids = ancestry_samples[i]
        }
        call Vcf2Counts as ancestry_all_counts {
            input:
                vcf_gz = ancestry_all.out_vcf_gz,
                tbi = ancestry_all.out_tbi
        }
        call Counts2Plot as ancestry_all_plot {
            input:
                gt_counts = ancestry_all_counts.gt_counts,
                out_file_name = ancestry_name[i]+"_all.png",
                plothw_r = plothw_r
        }
        
        # Frequent
        call FilterBySamples as ancestry_frequent {
            input:
                bcf = frequent.out_vcf_gz,
                csi = frequent.out_tbi,
                sample_ids = ancestry_samples[i]
        }
        call Vcf2Counts as ancestry_frequent_counts {
            input:
                vcf_gz = ancestry_frequent.out_vcf_gz,
                tbi = ancestry_frequent.out_tbi
        }
        call Counts2Plot as ancestry_frequent_plot {
            input:
                gt_counts = ancestry_frequent_counts.gt_counts,
                out_file_name = ancestry_name[i]+"_frequent.png",
                plothw_r = plothw_r
        }
    }
    
    output {
    }
}


# Remark: the task accepts both a VCF.GZ and a .BCF in input.
#
task SelectBiallelic {
    input {
        File vcf_gz
        File vcf_tbi
        
        Int n_cpu = 2
        Int ram_size_gb = 128
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 5*ceil(size(vcf_gz,"GB"))

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --max-alleles 2 --output-type z ~{vcf_gz} > biallelic.vcf.gz
        ${TIME_COMMAND} tabix -f biallelic.vcf.gz
    >>>

    output {
        File out_vcf_gz = "biallelic.vcf.gz"
        File out_tbi = "biallelic.vcf.gz.tbi"
    }

    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# Keeps only records with >= or < a given number of samples in which they were
# discovered.
#
task FilterByAf {
    input {
        File bcf
        File csi
        Float min_af
        Int smaller_or_larger
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        smaller_or_larger: "0: <min_af, 1: >=min_af"
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 3*ceil(size(bcf,"GB"))

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        if [ ~{smaller_or_larger} -eq 0 ]; then
            OPERATOR='<'
        else
            OPERATOR='>='
        fi
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include 'AF'${OPERATOR}~{min_af} --output-type z ~{bcf} > out.vcf.gz
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>

    output {
        File out_vcf_gz = "out.vcf.gz"
        File out_tbi = "out.vcf.gz.tbi"
    }

    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# The output of this program contains only records whose GT has two alleles and
# that are ALT in some sample.
#
task Vcf2Counts {
    input {
        File vcf_gz
        File tbi
        File? PlotHw_java
        
        Int n_cpu = 1
        Int ram_size_gb = 8
    }
    parameter_meta {
        PlotHw_java: "Custom Java program to use."
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 2*ceil(size(vcf_gz,"GB"))

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))

        
        if ~{defined(PlotHw_java)}
        then
            mv ~{PlotHw_java} ./PlotHwFast.java
            javac PlotHwFast.java
            CP_STRING=" "
        else
            CP_STRING="-cp ~{docker_dir}"
        fi
        ${TIME_COMMAND} java ${CP_STRING} -Xmx${EFFECTIVE_RAM_GB}G PlotHwFast ~{vcf_gz} gt_counts.csv
    >>>

    output {
        File gt_counts = "gt_counts.csv"
    }

    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


#
task Counts2Plot {
    input {
        File gt_counts
        String out_file_name
        File? plothw_r
    }
    parameter_meta {
        plothw_r: "Custom R script to use for plotting."
    }
    
    Int disk_size_gb = 10*ceil(size(gt_counts,"GB"))

    command <<<
        set -euxo pipefail

        if ~{defined(plothw_r)}
        then
            Rscript ~{plothw_r} ~{gt_counts} ~{out_file_name}.png
        else
            Rscript /hwe/PlotHW.r ~{gt_counts} ~{out_file_name}.png
        fi
    >>>

    output {
        File out_image = out_file_name + ".png"
    }

    runtime {
        docker: "fcunial/hapestry:hwe"
        cpu: 1
        memory: "4G"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}


# Remark: only samples that are present both in `sample_ids` and `vcf_gz` are
# kept, and only records that are ALT in some selected sample are kept.
#
task FilterBySamples {
    input {
        File bcf
        File csi
        File sample_ids
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        sample_ids: "One sample per line. Not necessarily sorted."
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 3*ceil(size(bcf,"GB"))

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        cut -f 1 ~{sample_ids} | sort > desired_samples.txt
        bcftools view --header-only ~{bcf} | tail -n 1 | tr '\t' '\n' | tail -n +10 | sort > present_samples.txt
        comm -1 -2 desired_samples.txt present_samples.txt > selected_samples.txt
        date
        bcftools view --threads ${N_THREADS} --samples-file selected_samples.txt ~{bcf} | bcftools filter --include 'COUNT(GT="alt")>0' --output-type z > out.vcf.gz
        date
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>

    output {
        File out_vcf_gz = "out.vcf.gz"
        File out_tbi = "out.vcf.gz.tbi"
        File selected_samples = "selected_samples.txt"
    }

    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
