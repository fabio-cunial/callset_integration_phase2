version 1.0

import "PlotHweImpl.wdl" as impl


#
workflow PlotHweFocused {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        Array[File] sample_ids
        
        Int max_distance_bp = 10
        
        File? plothw_r
    }
    
    call impl.FilterByLengthAndType as all_50 {
        input:
            vcf_gz = intersample_vcf_gz,
            vcf_tbi = intersample_tbi,
            min_sv_length = 50,
            sv_type = 0
    }
    
    # Filtering by ancestry
    scatter (samples in sample_ids) {
        call impl.FilterBySamples as ancestry {
            input:
                vcf_gz = all_50.out_vcf_gz,
                vcf_tbi = all_50.out_tbi,
                sample_ids = samples
        }
        call impl.Vcf2Counts as counts {
            input:
                vcf_gz = ancestry.out_vcf_gz,
                vcf_tbi = ancestry.out_tbi
        }
        call impl.Counts2Plot as plot {
            input:
                gt_counts = counts.gt_counts,
                out_file_name = "ancestry.png",
                plothw_r = plothw_r
        }
    }
    
    # Filtering out multiallelics
    call impl.SelectBiallelic as biallelic_50 {
        input:
            vcf_gz = all_50.out_vcf_gz,
            vcf_tbi = all_50.out_tbi,
            max_distance_bp = 10
    }
    call impl.Vcf2Counts as counts_biallelic {
        input:
            vcf_gz = biallelic_50.out_vcf_gz,
            vcf_tbi = biallelic_50.out_tbi
    }
    call impl.Counts2Plot as plot_biallelic {
        input:
            gt_counts = counts_biallelic.gt_counts,
            out_file_name = "ancestry.png",
            plothw_r = plothw_r
    }
    
    output {
    }
}
