version 1.0

#
workflow PlotHweImpl {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        File tandem_track_bed
        
        Int min_allele_count = 2
        Int max_distance_bp = 10
        
        File? plothw_r
    }
    
    # All
    call Vcf2Counts as counts {
        input:
            vcf_gz = intersample_vcf_gz,
            vcf_tbi = intersample_tbi
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
    call Counts2Plot as del_plot {
        input:
            gt_counts = del_counts.gt_counts,
            out_file_name = "del.png",
            plothw_r = plothw_r
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
    call Counts2Plot as ins_plot {
        input:
            gt_counts = ins_counts.gt_counts,
            out_file_name = "ins.png",
            plothw_r = plothw_r
    }
    
    # TRs
    call SelectTRs as trs {
        input:
            vcf_gz = intersample_vcf_gz,
            vcf_tbi = intersample_tbi,
            tandem_track_bed = tandem_track_bed,
            mode = 1
    }
    # Not TRs
    call SelectTRs as not_trs {
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
    call Counts2Plot as trs_plot {
        input:
            gt_counts = trs_counts.gt_counts,
            out_file_name = "trs.png",
            plothw_r = plothw_r
    }
    call Vcf2Counts as not_trs_counts {
        input:
            vcf_gz = not_trs.out_vcf_gz,
            vcf_tbi = not_trs.out_tbi
    }
    call Counts2Plot as not_trs_plot {
        input:
            gt_counts = not_trs_counts.gt_counts,
            out_file_name = "not_trs.png",
            plothw_r = plothw_r
    }
    
    # Frequent
    call FilterByAc as ac {
        input:
            vcf_gz = intersample_vcf_gz,
            vcf_tbi = intersample_tbi,
            min_count = 2
    }
    call Vcf2Counts as ac_counts {
        input:
            vcf_gz = ac.out_vcf_gz,
            vcf_tbi = ac.out_tbi
    }
    call Counts2Plot as ac_plot {
        input:
            gt_counts = ac_counts.gt_counts,
            out_file_name = "ac.png",
            plothw_r = plothw_r
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
    Int disk_size_gb = 3*ceil(size(vcf_gz,"GB"))

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
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --regions-file list.bed --include "SVLEN>=~{min_sv_length} || SVLEN<=-~{min_sv_length}" ~{vcf_gz} --output-type z > filtered.vcf.gz
        elif [ ~{sv_type} -eq 1 ]; then
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --regions-file list.bed --include "SVTYPE==\"DEL\" && (SVLEN>=~{min_sv_length} || SVLEN<=-~{min_sv_length})" ~{vcf_gz} --output-type z > filtered.vcf.gz
        elif [ ~{sv_type} -eq 2 ]; then
            ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --regions-file list.bed --include "SVTYPE==\"INS\" && (SVLEN>=~{min_sv_length} || SVLEN<=-~{min_sv_length})" ~{vcf_gz} --output-type z > filtered.vcf.gz
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

        
        source activate truvari5
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
        File vcf_gz
        File vcf_tbi
        Int min_count = 2
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 3*ceil(size(vcf_gz,"GB"))

    command <<<
        set -euxo pipefail
        
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
        File vcf_gz
        File vcf_tbi
        File tandem_track_bed
        Int mode
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        mode: "Keep only calls that: 0=do not overlap with the track; 1=overlap with the track."
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 10*ceil(size(vcf_gz,"GB"))

    command <<<
        set -euxo pipefail
        
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


#
task Vcf2Counts {
    input {
        File vcf_gz
        File vcf_tbi
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
        fi
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB}G PlotHwFast ~{vcf_gz} gt_counts.csv
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


#
task FilterBySamples {
    input {
        File vcf_gz
        File vcf_tbi
        File sample_ids
        
        Int n_cpu = 8
        Int ram_size_gb = 16
    }
    parameter_meta {
        sample_ids: "One sample per line. Not necessarily sorted."
    }
    
    String docker_dir = "/callset_integration"
    Int disk_size_gb = 3*ceil(size(vcf_gz,"GB"))

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        cut -f 1 ~{sample_ids} | sort > desired_samples.txt
        bcftools view --header-only ~{vcf_gz} | tail -n 1 | tr '\t' '\n' | tail -n +10 | sort > present_samples.txt
        comm -1 -2 desired_samples.txt present_samples.txt > selected_samples.txt
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file selected_samples.txt --output-type z ~{vcf_gz} > out.vcf.gz
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
