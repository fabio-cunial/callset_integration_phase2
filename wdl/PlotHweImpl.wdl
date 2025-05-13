version 1.0

#
workflow PlotHweImpl {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        File tandem_track_bed
        
        Int min_allele_count = 2
        Int max_distance_bp = 10
    }
    
    # All
    call Vcf2Counts as counts {
        input:
            vcf_gz = intersample_vcf_gz,
            vcf_tbi = intersample_tbi
    }
    call Counts2Plot as pdf {
        input:
            gt_counts = counts.gt_counts
    }
    
    # DEL
    call FilterByLengthAndType as del {
        input:
            vcf_gz = intersample_vcf_gz,
            vcf_tbi = intersample_tbi,
            min_sv_length = 1,
            sv_type = 1
    }
    call Vcf2Counts as del_counts {
        input:
            vcf_gz = del.out_vcf_gz,
            vcf_tbi = del.out_tbi
    }
    call Counts2Plot as del_pdf {
        input:
            gt_counts = del_counts.gt_counts
    }
    
    # INS
    call FilterByLengthAndType as ins {
        input:
            vcf_gz = intersample_vcf_gz,
            vcf_tbi = intersample_tbi,
            min_sv_length = 1,
            sv_type = 2
    }
    call Vcf2Counts as ins_counts {
        input:
            vcf_gz = ins.out_vcf_gz,
            vcf_tbi = ins.out_tbi
    }
    call Counts2Plot as ins_pdf {
        input:
            gt_counts = ins_counts.gt_counts
    }
    
    # Not in TRs
    call SelectTRs as trs {
        input:
            vcf_gz = intersample_vcf_gz,
            vcf_tbi = intersample_tbi,
            tandem_track_bed = tandem_track_bed,
            mode = 0
    }
    call Vcf2Counts as trs_counts {
        input:
            vcf_gz = trs.out_vcf_gz,
            vcf_tbi = trs.out_tbi
    }
    call Counts2Plot as trs_pdf {
        input:
            gt_counts = trs_counts.gt_counts
    }
    
    # Frequent
    call FilterByAc {
        input:
            vcf_gz = intersample_vcf_gz,
            vcf_tbi = intersample_tbi,
            min_count = 2
    }
    
    output {
    }
}


# Remark: the task keeps only autosomal calls.
#
task FilterByLengthAndType {
    input {
        File vcf_gz
        File vcf_tbi
        Int min_sv_length
        Int sv_type
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        sv_type: "0=all, 1=del, 2=ins."
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 3*ceil(size(vcf_gz,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        if [ ~{sv_type} -eq 0 ]; then
            INCLUDE_PREFIX=" "
            INCLUDE_STR=" "
        elif [ ~{sv_type} -eq 1 ]; then
            INCLUDE_PREFIX="--include"
            INCLUDE_STR="SVTYPE==\"DEL\" && (SVLEN>=~{min_sv_length} || SVLEN<=-~{min_sv_length})"
        elif [ ~{sv_type} -eq 2 ]; then
            INCLUDE_PREFIX="--include"
            INCLUDE_STR="SVTYPE==\"INS\" && (SVLEN>=~{min_sv_length} || SVLEN<=-~{min_sv_length})"
        fi
        for CHR in $(seq 1 22); do
            echo "chr${CHR}\t0\t3000000000" >> list.bed
        done
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --regions-file list.txt ${INCLUDE_PREFIX} ${INCLUDE_STR} ~{vcf_gz} --output-type z > filtered.vcf.gz
        ${TIME_COMMAND} tabix -f filtered.vcf.gz
    >>>

    output {
        File out_vcf_gz = work_dir + "/filtered.vcf.gz"
        File out_tbi = work_dir + "/filtered.vcf.gz.tbi"
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
        Int ram_size_gb = 16
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 5*ceil(size(vcf_gz,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        
        source activate truvari5
        date
        truvari anno numneigh --sizemin 1 --refdist ~{max_distance_bp} ~{vcf_gz} | bcftools view --include 'INFO/NumNeighbors == 0' | bgzip --compression-level 2 > biallelic.vcf.gz
        date
        ${TIME_COMMAND} tabix -f biallelic.vcf.gz
        date
    >>>

    output {
        File out_vcf_gz = work_dir + "/biallelic.vcf.gz"
        File out_tbi = work_dir + "/biallelic.vcf.gz.tbi"
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
        File vcf_gz
        File vcf_tbi
        Int min_count = 2
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 3*ceil(size(vcf_gz,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "AC>=~{min_count}" --output-type z ~{vcf_gz} > out.vcf.gz
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>

    output {
        File out_vcf_gz = work_dir + "/out.vcf.gz"
        File out_tbi = work_dir + "/out.vcf.gz.tbi"
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
        File vcf_gz
        File vcf_tbi
        File tandem_track_bed
        Int mode
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        mode: "0=outside the track; 1=inside the track."
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 10*ceil(size(vcf_gz,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        bcftools view --header-only ~{vcf_gz} > out.vcf
        if [ ~{mode} -eq 1 ]; then
            INTERSECTION_MODE="-u"
        elif [ ~{mode} -eq 0 ]; then
            INTERSECTION_MODE="-v"
        fi
        ${TIME_COMMAND} bedtools intersect -a ~{vcf_gz} -b ~{tandem_track_bed} ${INTERSECTION_MODE} >> out.vcf
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} --compress-level 2 out.vcf
        ${TIME_COMMAND} tabix -f out.vcf.gz
    >>>

    output {
        File out_vcf_gz = work_dir + "/out.vcf.gz"
        File out_tbi = work_dir + "/out.vcf.gz.tbi"
    }

    runtime {
        docker: "fcunial/callset_integration_phase2_workpackages"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}


#
task Vcf2Counts {
    input {
        File vcf_gz
        File vcf_tbi
        
        Int n_cpu = 1
        Int ram_size_gb = 8
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 2*ceil(size(vcf_gz,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 1 ))

        
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB}G PlotHwFast ~{vcf_gz} gt_counts.csv
    >>>

    output {
        File gt_counts = work_dir + "/gt_counts.csv"
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
    }
    
    Int disk_size_gb = 10*ceil(size(gt_counts,"GB"))

    command <<<
        set -euxo pipefail

        Rscript /hwe/PlotHW.r ~{gt_counts} hwe_plot.pdf
    >>>

    output {
        File out_pdf = "hwe_plot.pdf"
    }

    runtime {
        docker: "fcunial/hapestry:hwe"
        cpu: 1
        memory: "4G"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}
