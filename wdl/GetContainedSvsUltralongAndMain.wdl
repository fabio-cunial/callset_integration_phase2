version 1.0


#
workflow GetContainedSvsUltralongAndMain {
    input {
        File ultralong_bcf
        File ultralong_csi
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
            ultralong_csi = ultralong_csi,
            main_bcf = main_bcf,
            main_csi = main_csi,
            min_calls = min_calls,
            min_sv_length = min_sv_length,
            docker_image = docker_image
    }
    
    output {
        File out_bed = Impl.out_bed
    }
}


task Impl {
    input {
        File ultralong_bcf
        File ultralong_csi
        File main_bcf
        File main_csi

        Int min_calls
        Int min_sv_length
        
        String docker_image
        Int n_cpu = 4
        Int ram_size_gb = 8
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil( size(ultralong_bcf,"GB") ) + 10*ceil( size(main_bcf,"GB") )
    String docker_dir = "/callset_integration"
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))




        # ----------------------- Steps of the pipeline ------------------------

cat << 'END' > analyze_interval.sh
#!/bin/bash
set -euo pipefail

CONTAINER_SVTYPE=$1
ID=$2
CONTAINER_ROW=$3

CHROM=$(echo ${CONTAINER_ROW} | cut -d , -f 1)
START=$(echo ${CONTAINER_ROW} | cut -d , -f 2)
END=$(echo ${CONTAINER_ROW} | cut -d , -f 3)

echo -e "${CHROM}\t${START}\t${END}" > ${ID}.bed
bcftools view --regions-file ${ID}.bed --regions-overlap variant --output-type v ~{main_bcf} --output ${ID}.vcf
rm -f ${ID}.bed
N_CONTAINED_RECORDS=$( grep -c -v '^#' ${ID}.vcf || true )
if [ ${N_CONTAINED_RECORDS} -ge ~{min_calls} ]; then
    java -cp ~{docker_dir} GetContainedSvsUltralongAndMain ${ID}.vcf ${N_CONTAINED_RECORDS} ${CONTAINER_SVTYPE} ~{min_calls} ~{min_sv_length} ${ID}_counts.bed ${CONTAINER_ROW}
fi
rm -f ${ID}.vcf
END
        chmod +x analyze_interval.sh


        function GetContainedSvs() {
            local SVTYPE=$1

            date 1>&2
            if [ ${SVTYPE} = "INSDUP" ]; then
                (bcftools view --header-only ~{ultralong_bcf} ; bcftools view --no-header --include 'SVTYPE="INS"' ~{ultralong_bcf} | { grep INSDUP || true; }) | bcftools query --format '%CHROM,%INSDUP_POS,%INSDUP_SVLEN,[%GT,]\n' | awk 'BEGIN { FS=","; OFS=","; i=0; } { $3=$2+$3; printf("%d %s\n",i++,$0); }' > ${SVTYPE}_intervals.wsv
            else
                bcftools query --format '%CHROM,%POS,%SVLEN,[%GT,]\n' --include 'SVTYPE="'${SVTYPE}'"' ~{ultralong_bcf} | awk 'BEGIN { FS=","; OFS=","; i=0; } { $3=$2+$3; printf("%d %s\n",i++,$0); }' > ${SVTYPE}_intervals.wsv
            fi
            if [ ! -s ${SVTYPE}_intervals.wsv ]; then
                touch ${SVTYPE}_out.bed
                return
            fi
            date 1>&2
            ${TIME_COMMAND} xargs --arg-file=${SVTYPE}_intervals.wsv --max-lines=1 --max-procs=${N_THREADS} ./analyze_interval.sh ${SVTYPE}
            rm -f ${SVTYPE}_intervals.wsv
            find . -maxdepth 1 -name '*_counts.bed' -exec cat {} + > ${SVTYPE}_out.bed
            find . -maxdepth 1 -name '*_counts.bed' -delete
        }


        function GetContainerSvs() {
            # To be implemented...
            :
        }




        # ---------------------------- Main program ----------------------------

        # Intervals
        GetContainedSvs DEL
        GetContainedSvs INSDUP
        GetContainedSvs DUP
        GetContainedSvs INV
        
        # Points
        GetContainerSvs INS

        # Outputting
        cat *_out.bed > out.bed
    >>>
    
    output {
        File out_bed = "out.bed"
    }
    runtime {
        docker: docker_image
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
