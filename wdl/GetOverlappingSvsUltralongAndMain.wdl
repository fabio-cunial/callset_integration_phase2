version 1.0


#
workflow GetOverlappingSvsUltralongAndMain {
    input {
        File ultralong_bcf
        File main_bcf
        File main_csi

        Int min_calls
        Int min_sv_length

        String docker_image = "us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages"
    }
    parameter_meta {
    }
    
    call Impl {
        input:
            ultralong_bcf = ultralong_bcf,
            main_bcf = main_bcf,
            main_csi = main_csi,
            min_calls = min_calls,
            min_sv_length = min_sv_length,
            docker_image = docker_image
    }
    
    output {
        File out_intervals_bed = Impl.out_intervals_bed
        File out_ins_bed = Impl.out_ins_bed
    }
}


task Impl {
    input {
        File ultralong_bcf
        File main_bcf
        File main_csi

        Int min_calls
        Int min_sv_length
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
        Int disk_size_gb = 1000
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        RAM_PER_THREAD_MB=$(( (~{ram_size_gb} * 1024 - 500) / ${N_THREADS} ))




        # ----------------------- Steps of the pipeline ------------------------

cat << 'END' > process_chunk.sh
#!/bin/bash
set -euxo pipefail

CONTAINER_SVTYPE=$1
CHUNK_ID=$2
CHUNK_FILE=$3
RAM_PER_THREAD_MB=$4

# Creating contained VCFs for every container call
N_CALLS=$(wc -l < ${CHUNK_FILE})
CALL_ID="0"; QUANTUM="100"
while IFS=, read -u 3 CHROM START END REST; do
    if [ $(( ${CALL_ID} % ${QUANTUM} )) -eq 0 ]; then
        df -h 1>&2
    fi
    echo -e "${CHROM}\t${START}\t${END}" > ${CHUNK_ID}_${CALL_ID}.bed
    cat header.tsv > ${CHUNK_ID}_${CALL_ID}.vcf
    bcftools view --no-header --regions-file ${CHUNK_ID}_${CALL_ID}.bed --regions-overlap variant --output-type v ~{main_bcf} >> ${CHUNK_ID}_${CALL_ID}.vcf
    rm -f ${CHUNK_ID}_${CALL_ID}.bed
    CALL_ID=$(( ${CALL_ID} + 1 ))
done 3< ${CHUNK_FILE}

# Processing all container calls and contained VCFs
/usr/bin/time --verbose java -cp ~{docker_dir} -Xmx${RAM_PER_THREAD_MB}M GetOverlappingSvsUltralongAndMain ${CONTAINER_SVTYPE} ~{min_calls} ~{min_sv_length} ${CHUNK_ID} ${CHUNK_FILE} . ${CHUNK_ID}_counts.bed
rm -f ${CHUNK_ID}_*.vcf
END
        chmod +x process_chunk.sh


        # Processes every ultralong call of a given interval type, looking up 
        # overlapping calls in the main VCF.
        #
        function GetOverlappingSvs() {
            local SVTYPE=$1

            date 1>&2
            if [ ${SVTYPE} = "INSDUP" ]; then
                (bcftools view --header-only ~{ultralong_bcf} ; bcftools view --no-header --include 'SVTYPE="INS"' ~{ultralong_bcf} | { grep INSDUP || true; }) | bcftools query --format '%CHROM,%INSDUP_POS,%INSDUP_SVLEN,[%GT,]\n' | awk 'BEGIN { FS=","; OFS=","; } { $3 = $2 + ($3 < 0 ? -$3 : $3); printf("%s\n",$0); }' > ${SVTYPE}_intervals.wsv
            elif [ ${SVTYPE} = "INS" ]; then
                (bcftools view --header-only ~{ultralong_bcf} ; bcftools view --no-header --include 'SVTYPE="INS"' ~{ultralong_bcf} | { grep -v INSDUP || true; }) | bcftools query --format '%CHROM,%POS,%POS,[%GT,]\n' | awk 'BEGIN { FS=","; OFS=","; } { $2 = $2 - 5; $3 = $3 + 5; printf("%s\n",$0); }' > ${SVTYPE}_intervals.wsv
            else
                bcftools query --format '%CHROM,%POS,%SVLEN,[%GT,]\n' --include 'SVTYPE="'${SVTYPE}'"' ~{ultralong_bcf} | awk 'BEGIN { FS=","; OFS=","; } { $3 = $2 + ($3 < 0 ? -$3 : $3); printf("%s\n",$0); }' > ${SVTYPE}_intervals.wsv
            fi
            date 1>&2
            if [ ! -s ${SVTYPE}_intervals.wsv ]; then
                touch ${SVTYPE}_out.bed
                return
            fi
            N_INTERVALS=$(wc -l < ${SVTYPE}_intervals.wsv)
            N_INTERVALS_PER_THREAD=$(( ( ${N_INTERVALS} + ${N_THREADS} - 1 ) / ${N_THREADS} ))
            split -d -a 2 -l ${N_INTERVALS_PER_THREAD} ${SVTYPE}_intervals.wsv ${SVTYPE}_chunk_
            rm -f ${SVTYPE}_intervals.wsv
            local PIDS=()
            for CHUNK_FILE in ${SVTYPE}_chunk_*; do
                ID=$(basename ${CHUNK_FILE} | cut -d _ -f 3)
                ${TIME_COMMAND} ./process_chunk.sh ${SVTYPE} ${ID} ${CHUNK_FILE} ${RAM_PER_THREAD_MB} &
                PIDS+=($!)
            done
            for P in "${PIDS[@]}"; do wait ${P}; done
            rm -f ${SVTYPE}_chunk_*
            cat *_counts.bed > ${SVTYPE}_out.bed
            rm -f *_counts.bed
        }




        # ---------------------------- Main program ----------------------------
        
        bcftools view --header-only ~{main_bcf} | tail -n 1 > header.tsv

        GetOverlappingSvs DEL
        GetOverlappingSvs INSDUP
        GetOverlappingSvs DUP
        GetOverlappingSvs INV
        GetOverlappingSvs INS

        # Sorting by decreasing density
        cat DEL_out.bed INSDUP_out.bed DUP_out.bed INV_out.bed | sort -k6,6gr > out_intervals.bed
        sort -k5,5nr INS_out.bed > out_ins.bed
    >>>
    
    output {
        File out_intervals_bed = "out_intervals.bed"
        File out_ins_bed = "out_ins.bed"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
