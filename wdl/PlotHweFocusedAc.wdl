version 1.0

import "PlotHweImpl.wdl" as impl


#
workflow PlotHweFocusedAc {
    input {
        File intersample_vcf_gz
        File intersample_tbi
        Array[File] sample_ids
        Array[String] ancestry_ids
        
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
        # Ancestry, AC.
        scatter (j in range(10)) {
            call impl.FilterByAc as ac {
                input:
                    vcf_gz = ancestry.out_vcf_gz,
                    vcf_tbi = ancestry.out_tbi,
                    min_count = j
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
                    out_file_name = ancestry_ids[i] + "_ac" + j,
                    plothw_r = plothw_r
            }
        }
    }
    
    output {
    }
}
