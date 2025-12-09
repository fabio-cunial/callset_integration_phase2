version 1.0

#
workflow SV_Integration_PlotHwe {
    input {
        File intersample_bcf
        File intersample_csi
        File tandem_track_bed
        
        Int min_allele_count = 2
        Int max_distance_bp = 10
        
        File? plothw_r
    }
    
    # All
    call Bcf2Vcf {
        input:
            bcf = intersample_bcf,
            csi = intersample_csi
    }
    call Vcf2Counts as counts {
        input:
            vcf_gz = Bcf2Vcf.out_vcf_gz,
            tbi = Bcf2Vcf.out_tbi
    }
    call Counts2Plot as plot {
        input:
            gt_counts = counts.gt_counts,
            out_file_name = "all.png",
            plothw_r = plothw_r
    }
    
    # DEL
    call FilterByLengthAndType as del {
        input:
            bcf = intersample_bcf,
            csi = intersample_csi,
            min_sv_length = 1,
            sv_type = 1
    }
    call Vcf2Counts as del_counts {
        input:
            vcf_gz = del.out_vcf_gz,
            tbi = del.out_tbi
    }
    call Counts2Plot as del_plot {
        input:
            gt_counts = del_counts.gt_counts,
            out_file_name = "del.png",
            plothw_r = plothw_r
    }
    
    # INS
    call FilterByLengthAndType as ins {
        input:
            bcf = intersample_bcf,
            csi = intersample_csi,
            min_sv_length = 1,
            sv_type = 2
    }
    call Vcf2Counts as ins_counts {
        input:
            vcf_gz = ins.out_vcf_gz,
            tbi = ins.out_tbi
    }
    call Counts2Plot as ins_plot {
        input:
            gt_counts = ins_counts.gt_counts,
            out_file_name = "ins.png",
            plothw_r = plothw_r
    }
    
    # TRs
    call SelectTRs as trs {
        input:
            bcf = intersample_bcf,
            csi = intersample_csi,
            tandem_track_bed = tandem_track_bed,
            mode = 1
    }
    call Vcf2Counts as trs_counts {
        input:
            vcf_gz = trs.out_vcf_gz,
            tbi = trs.out_tbi
    }
    call Counts2Plot as trs_plot {
        input:
            gt_counts = trs_counts.gt_counts,
            out_file_name = "trs.png",
            plothw_r = plothw_r
    }
    
    # Not TRs
    call SelectTRs as not_trs {
        input:
            bcf = intersample_bcf,
            csi = intersample_csi,
            tandem_track_bed = tandem_track_bed,
            mode = 0
    }
    call Vcf2Counts as not_trs_counts {
        input:
            vcf_gz = not_trs.out_vcf_gz,
            tbi = not_trs.out_tbi
    }
    call Counts2Plot as not_trs_plot {
        input:
            gt_counts = not_trs_counts.gt_counts,
            out_file_name = "not_trs.png",
            plothw_r = plothw_r
    }
    
    # Frequently discovered
    call FilterByNDiscoverySamplesGeq as frequent {
        input:
            bcf = intersample_bcf,
            csi = intersample_csi,
            min_count = 2
    }
    call Vcf2Counts as frequent_counts {
        input:
            vcf_gz = frequent.out_vcf_gz,
            tbi = frequent.out_tbi
    }
    call Counts2Plot as frequent_plot {
        input:
            gt_counts = frequent_counts.gt_counts,
            out_file_name = "ac.png",
            plothw_r = plothw_r
    }
    
    # Infrequently discovered
    call FilterByNDiscoverySamplesLt as infrequent {
        input:
            bcf = intersample_bcf,
            csi = intersample_csi,
            min_count = 2
    }
    call Vcf2Counts as infrequent_counts {
        input:
            vcf_gz = infrequent.out_vcf_gz,
            tbi = infrequent.out_tbi
    }
    call Counts2Plot as infrequent_plot {
        input:
            gt_counts = infrequent_counts.gt_counts,
            out_file_name = "ac.png",
            plothw_r = plothw_r
    }
    
    output {
    }
}


#
task Bcf2Vcf {
    input {
        File bcf
        File csi
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 3*ceil(size(bcf,"GB"))

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --output-type z ~{bcf} > out.vcf.gz
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


# Remark: the task keeps only autosomal calls.
#
task FilterByLengthAndType {
    input {
        File bcf
        File csi
        Int min_sv_length
        Int sv_type
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        sv_type: "0=all, 1=del, 2=ins."
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 3*ceil(size(bcf,"GB"))

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        for CHR in $(seq 1 22); do
            echo -e "chr${CHR}\t0\t3000000000" >> list.bed
        done
        if [ ~{sv_type} -eq 0 ]; then
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --regions-file list.bed --include "ABS(SVLEN)>=~{min_sv_length}" --output-type z ~{bcf} > filtered.vcf.gz
        elif [ ~{sv_type} -eq 1 ]; then
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --regions-file list.bed --include "SVTYPE==\"DEL\" && ABS(SVLEN)>=~{min_sv_length}" --output-type z ~{bcf} > filtered.vcf.gz
        elif [ ~{sv_type} -eq 2 ]; then
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --regions-file list.bed --include "SVTYPE==\"INS\" && ABS(SVLEN)>=~{min_sv_length}" --output-type z ~{bcf} > filtered.vcf.gz
        fi
        ${TIME_COMMAND} tabix -f filtered.vcf.gz
    >>>

    output {
        File out_vcf_gz = "filtered.vcf.gz"
        File out_tbi = "filtered.vcf.gz.tbi"
    }

    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


# A biallelic call is assumed to be an isolated call with no neighbor within 
# `max_distance_bp`.
#
task SelectBiallelic {
    input {
        File vcf_gz
        File vcf_tbi
        Int max_distance_bp = 10
        
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

        
        date
        truvari anno numneigh --sizemin 1 --refdist ~{max_distance_bp} ~{vcf_gz} | bcftools view --include 'INFO/NumNeighbors == 0' | bgzip -@ ${N_THREADS} --compress-level 2 > biallelic.vcf.gz
        date
        ${TIME_COMMAND} tabix -f biallelic.vcf.gz
        date
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


# Keeps only ALT alleles with at least a given count.
#
task FilterByAc {
    input {
        File bcf
        File csi
        Int min_count = 2
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
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
        
        
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "AC>=~{min_count}" --output-type z ~{bcf} > out.vcf.gz
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


# Keeps only records with >= a given number of samples in which they were
# discovered.
#
task FilterByNDiscoverySamplesGeq {
    input {
        File bcf
        File csi
        Int min_count = 1268
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
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
        
        
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "N_DISCOVERY_SAMPLES>=~{min_count}" --output-type z ~{bcf} > out.vcf.gz
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


# Keeps only records with < a given number of samples in which they were
# discovered.
#
task FilterByNDiscoverySamplesLt {
    input {
        File bcf
        File csi
        Int min_count = 1268
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
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
        
        
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "N_DISCOVERY_SAMPLES<~{min_count}" --output-type z ~{bcf} > out.vcf.gz
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


# Any overlap with the track is considered, even by a single bp.
#
task SelectTRs {
    input {
        File bcf
        File csi
        File tandem_track_bed
        Int mode
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        mode: "Keep only calls that: 0=do not overlap with the track; 1=overlap with the track."
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 10*ceil(size(bcf,"GB"))

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        bcftools view --header-only ~{bcf} > out.vcf
        if [ ~{mode} -eq 1 ]; then
            INTERSECTION_MODE="-u"
        elif [ ~{mode} -eq 0 ]; then
            INTERSECTION_MODE="-v"
        fi
        ${TIME_COMMAND} bedtools intersect -a ~{bcf} -b ~{tandem_track_bed} ${INTERSECTION_MODE} >> out.vcf
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} --compress-level 2 out.vcf
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
        
        cut -f 1 ~{sample_ids} | sort -V > desired_samples.txt
        bcftools view --header-only ~{bcf} | tail -n 1 | tr '\t' '\n' | tail -n +10 | sort -V > present_samples.txt
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
