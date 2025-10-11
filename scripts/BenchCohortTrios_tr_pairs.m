SOURCE_DIR='.';
FONT_SIZE=14;
COLORMAP='cubehelix';

ML_IDS=[1 2 4 8 16 32 64 128];
NR_IDS=[1 2 4 8 16 32 64];
IL_IDS=[2 4 6 8 10 12 14 16 18 20 40 80 160 320];

TOTAL_RECORDS_07_50BP=1018505;  % All records, of any length: 2816177




figure(1);

subplot(2,3,1);
A=load(sprintf('%s/ml_il_mendelian.csv',SOURCE_DIR));
imagesc(A.'); colormap(COLORMAP); colorbar; axis square;
xlabel('motif length'); xticks([1:length(ML_IDS)]); xticklabels(ML_IDS);
ylabel('interval length'); yticks([1:length(IL_IDS)]); yticklabels(IL_IDS);
title('Mendelian error (avg over 5 trios)');

subplot(2,3,2);
A=load(sprintf('%s/nr_il_mendelian.csv',SOURCE_DIR));
imagesc(A.'); colormap(COLORMAP); colorbar; axis square;
xlabel('n repeats in ref'); xticks([1:length(NR_IDS)]); xticklabels(NR_IDS);
ylabel('interval length'); yticks([1:length(IL_IDS)]); yticklabels(IL_IDS);
title('Mendelian error (avg over 5 trios)');

subplot(2,3,3);
A=load(sprintf('%s/ml_nr_mendelian.csv',SOURCE_DIR));
imagesc(A); colormap(COLORMAP); colorbar; axis square;
ylabel('motif length'); yticks([1:length(ML_IDS)]); yticklabels(ML_IDS);
xlabel('n repeats in ref'); xticks([1:length(NR_IDS)]); xticklabels(NR_IDS);
title('Mendelian error (avg over 5 trios)');








figure(2);

subplot(2,3,4);
A=load(sprintf('%s/ml_il_counts.csv',SOURCE_DIR));
imagesc((A.')./TOTAL_RECORDS_07_50BP); colormap('hot'); colorbar; axis square;
xlabel('motif length'); xticks([1:length(ML_IDS)]); xticklabels(ML_IDS);
ylabel('interval length'); yticks([1:length(IL_IDS)]); yticklabels(IL_IDS);
title('Fraction of all calls \geq 50bp');

subplot(2,3,5);
A=load(sprintf('%s/nr_il_counts.csv',SOURCE_DIR));
imagesc((A.')./TOTAL_RECORDS_07_50BP); colormap('hot'); colorbar; axis square;
xlabel('n repeats in ref'); xticks([1:length(NR_IDS)]); xticklabels(NR_IDS);
ylabel('interval length'); yticks([1:length(IL_IDS)]); yticklabels(IL_IDS);
title('Fraction of all calls \geq 50bp');

subplot(2,3,6);
A=load(sprintf('%s/ml_nr_counts.csv',SOURCE_DIR));
imagesc(A./TOTAL_RECORDS_07_50BP); colormap('hot'); colorbar; axis square;
ylabel('motif length'); yticks([1:length(ML_IDS)]); yticklabels(ML_IDS);
xlabel('n repeats in ref'); xticks([1:length(NR_IDS)]); xticklabels(NR_IDS);
title('Fraction of all calls \geq 50bp');




