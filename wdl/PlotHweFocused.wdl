version 1.0

import "PlotHweImpl.wdl" as impl


#
workflow PlotHweFocused {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        Array[File] sample_ids
        Array[String] ancestry_ids
        
        Int max_distance_bp = 10
        Int min_count = 2
        
        File? PlotHw_java
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
    scatter (i in range(length(sample_ids))) {
        # Ancestry, all.
        call impl.FilterBySamples as ancestry {
            input:
                vcf_gz = all_50.out_vcf_gz,
                vcf_tbi = all_50.out_tbi,
                sample_ids = sample_ids[i]
        }
        call impl.Vcf2Counts as counts {
            input:
                vcf_gz = ancestry.out_vcf_gz,
                vcf_tbi = ancestry.out_tbi,
                PlotHw_java = PlotHw_java
        }
        call impl.Counts2Plot as plot {
            input:
                gt_counts = counts.gt_counts,
                out_file_name = ancestry_ids[i],
                plothw_r = plothw_r
        }
        # Ancestry, biallelic.
        call impl.SelectBiallelic as biallelic {
            input:
                vcf_gz = ancestry.out_vcf_gz,
                vcf_tbi = ancestry.out_tbi,
                max_distance_bp = max_distance_bp
        }
        call impl.Vcf2Counts as counts_biallelic {
            input:
                vcf_gz = biallelic.out_vcf_gz,
                vcf_tbi = biallelic.out_tbi,
                PlotHw_java = PlotHw_java
        }
        call impl.Counts2Plot as plot_biallelic {
            input:
                gt_counts = counts_biallelic.gt_counts,
                out_file_name = ancestry_ids[i] + "_biallelic",
                plothw_r = plothw_r
        }
        # Ancestry, AC>1.
        call impl.FilterByAc as ac {
            input:
                vcf_gz = ancestry.out_vcf_gz,
                vcf_tbi = ancestry.out_tbi,
                min_count = min_count
        }
        call impl.Vcf2Counts as counts_ac {
            input:
                vcf_gz = ac.out_vcf_gz,
                vcf_tbi = ac.out_tbi,
                PlotHw_java = PlotHw_java
        }
        call impl.Counts2Plot as plot_ac {
            input:
                gt_counts = counts_ac.gt_counts,
                out_file_name = ancestry_ids[i] + "_biallelic",
                plothw_r = plothw_r
        }
    }
    
    # All, biallelic.
    call impl.SelectBiallelic as biallelic_50 {
        input:
            vcf_gz = all_50.out_vcf_gz,
            vcf_tbi = all_50.out_tbi,
            max_distance_bp = max_distance_bp
    }
    call impl.Vcf2Counts as counts_biallelic_50 {
        input:
            vcf_gz = biallelic_50.out_vcf_gz,
            vcf_tbi = biallelic_50.out_tbi,
            PlotHw_java = PlotHw_java
    }
    call impl.Counts2Plot as plot_biallelic_50 {
        input:
            gt_counts = counts_biallelic_50.gt_counts,
            out_file_name = "biallelic",
            plothw_r = plothw_r
    }
    
    output {
    }
}
