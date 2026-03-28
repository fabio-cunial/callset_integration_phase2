! pip install scikit-allel

import numpy as np
import tqdm
import allel
import sklearn.metrics
import matplotlib.pyplot as plt
import glob

# Downloading VCFs (can be skipped if already downloaded)
! gsutil -m cp "gs://______________________________________________________/scratch/cunial_intersample_vcf/v3/ultralong_annotate/merged/scored/score.vcf.gz*" ./scored-vcfs/
scored_vcfs = glob.glob(f'scored-vcfs/score.vcf.gz')
	 
# Plotting logic
def load_callsets(scored_vcfs):
    callsets = []
    for i, scored_vcf in tqdm.tqdm(enumerate(scored_vcfs)):
        callset = allel.read_vcf(scored_vcf, fields=['variants/SVTYPE',
                                                     'variants/SVLEN',
                                                     'variants/SCORE',
                                                     'variants/CALIBRATION_SENSITIVITY',
                                                     'variants/training',
                                                     'variants/extracted',
                                                     'variants/CHROM'])
        callsets.append(callset)
    return callsets

def plot_roc(callsets, calibration_sensitivity_thresholds=[0.7,0.9]):
    # Masks
    # Remark: we trained on chr6..Y. We test here on chr1..5.
    test_mask_lambda = lambda callset: ((callset['variants/CHROM'] == 'chr1') | 
                                        (callset['variants/CHROM'] == 'chr2') | 
                                        (callset['variants/CHROM'] == 'chr3') | 
                                        (callset['variants/CHROM'] == 'chr4') | 
                                        (callset['variants/CHROM'] == 'chr5'))
    
    fig = plt.figure(figsize=(5, 5))
    for mask_title, mask_color, additional_mask_lambda in [
                                               ['all', 'blue',
                                                lambda callset: test_mask_lambda(callset) ]]:
        
        for i, callset in tqdm.tqdm(enumerate(callsets)):
            additional_mask = additional_mask_lambda(callset)
            
            # Score ROC
            metrics_mask = ~np.isnan(callset['variants/SCORE']) & additional_mask
            fpr, tpr, thresholds = sklearn.metrics.roc_curve(
                y_true=callset['variants/training'][metrics_mask], 
                y_score=callset['variants/SCORE'][metrics_mask])
            plt.plot(fpr, tpr, label=mask_title if i==0 else None, c=mask_color, alpha=0.25)
            
            # CAL_SENS thresholds
            metrics_mask = ~np.isnan(callset['variants/CALIBRATION_SENSITIVITY']) & additional_mask
            fpr, tpr, thresholds = sklearn.metrics.roc_curve(
                y_true=callset['variants/training'][metrics_mask], 
                y_score=1 - callset['variants/CALIBRATION_SENSITIVITY'][metrics_mask])
            fpr_selection = fpr[[np.argmax(tpr >= cst) for cst in calibration_sensitivity_thresholds]]
            tpr_selection = tpr[[np.argmax(tpr >= cst) for cst in calibration_sensitivity_thresholds]]
            plt.scatter(fpr_selection, tpr_selection, label=None, 
                        c=mask_color, marker='d', s=4)

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'Deletions >10kb, chr1..5, HPRC+HGSVC, GRCh38.')
    plt.legend()
    plt.show()
    
    
callsets = load_callsets(scored_vcfs)
plot_roc(callsets)

