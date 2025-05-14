version 1.0

import "PlotHweImpl.wdl" as impl


#
workflow PlotHwe {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        File tandem_track_bed
        
        Int min_allele_count = 2
        Int max_distance_bp = 10
        
        File? plothw_r
    }
    
    # Main categories
    call impl.FilterByLengthAndType as all {
        input:
            vcf_gz = intersample_vcf_gz,
            vcf_tbi = intersample_tbi,
            min_sv_length = 1,
            sv_type = 0
    }
    call impl.FilterByLengthAndType as all_50 {
        input:
            vcf_gz = all.out_vcf_gz,
            vcf_tbi = all.out_tbi,
            min_sv_length = 50,
            sv_type = 0
    }
    call impl.SelectBiallelic as biallelic {
        input:
            vcf_gz = all.out_vcf_gz,
            vcf_tbi = all.out_tbi,
            max_distance_bp = 10
    }
    call impl.FilterByLengthAndType as biallelic_50 {
        input:
            vcf_gz = biallelic.out_vcf_gz,
            vcf_tbi = biallelic.out_tbi,
            min_sv_length = 50,
            sv_type = 0
    }
    
    # Analyzing each main category
    call impl.PlotHweImpl as all_impl {
        input:
            intersample_vcf_gz = all.out_vcf_gz,
            intersample_tbi = all.out_tbi,
            tandem_track_bed = tandem_track_bed,
            min_allele_count = 2,
            max_distance_bp = 10,
            plothw_r = plothw_r
    }
    call impl.PlotHweImpl as all_50_impl {
        input:
            intersample_vcf_gz = all_50.out_vcf_gz,
            intersample_tbi = all_50.out_tbi,
            tandem_track_bed = tandem_track_bed,
            min_allele_count = 2,
            max_distance_bp = 10,
            plothw_r = plothw_r
    }
    call impl.PlotHweImpl as biallelic_impl {
        input:
            intersample_vcf_gz = biallelic.out_vcf_gz,
            intersample_tbi = biallelic.out_tbi,
            tandem_track_bed = tandem_track_bed,
            min_allele_count = 2,
            max_distance_bp = 10,
            plothw_r = plothw_r
    }
    call impl.PlotHweImpl as biallelic_50_impl {
        input:
            intersample_vcf_gz = biallelic_50.out_vcf_gz,
            intersample_tbi = biallelic_50.out_tbi,
            tandem_track_bed = tandem_track_bed,
            min_allele_count = 2,
            max_distance_bp = 10,
            plothw_r = plothw_r
    }
    
    output {
    }
}
