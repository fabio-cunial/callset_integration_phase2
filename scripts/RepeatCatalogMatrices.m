INPUT_DIR='/Users/fcunial/git/callset_integration_phase2/scripts';
FONT_SIZE=14;


intervalLength_motifLength=load(sprintf('%s/intervalLength_motifLength.csv',INPUT_DIR));
intervalLength_nRepeatsInRef=load(sprintf('%s/intervalLength_nRepeatsInRef.csv',INPUT_DIR));
motifLength_nRepeatsInRef=load(sprintf('%s/motifLength_nRepeatsInRef.csv',INPUT_DIR));

subplot(1,3,1);
imagesc(log10(intervalLength_motifLength)); xlabel('motif length'); ylabel('interval length'); axis square; set(gca,'fontsize',FONT_SIZE); colormap hot;

subplot(1,3,2);
imagesc(log10(intervalLength_nRepeatsInRef)); ylabel('interval length'); xlabel('n repeats in ref'); axis square; set(gca,'fontsize',FONT_SIZE); colormap hot;

subplot(1,3,3)
imagesc(log10(motifLength_nRepeatsInRef)); ylabel('motif length'); xlabel('n repeats in ref'); axis square; set(gca,'fontsize',FONT_SIZE); colormap hot;
