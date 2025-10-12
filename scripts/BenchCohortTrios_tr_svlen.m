SOURCE_DIR='.';
DELTA=0.4;
FONT_SIZE=10;

SV_LENGTHS=[50,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,1000000];
LABELS={"50","100","200","300","400","500","600","700","800","900","1k","2k","3k","4k","5k","6k","7k","8k","9k","10k","1m"};
N_LENGTHS=length(SV_LENGTHS);



figure(1);

subplot(2,2,1); hold on;
A=load(sprintf('%s/tr_svlen_del.csv',SOURCE_DIR));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    for j=[1:N_LENGTHS]
        X=ones(nrows,1).*j -DELTA/2 + rand(nrows,1).*DELTA;
        Y=A(:,2*j)./(A(:,2*j)+A(:,2*j-1));
        plot(X,Y,'.');
    endfor
endfor
xlabel('SVLEN'); xticks([1:N_LENGTHS]); xticklabels(LABELS); title('DEL inside TR / Mendelian error'); axis square; grid on; set(gca,'fontsize',FONT_SIZE);

subplot(2,2,2); hold on;
A=load(sprintf('%s/tr_svlen_ins.csv',SOURCE_DIR));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    for j=[1:N_LENGTHS]
        X=ones(nrows,1).*j -DELTA/2 + rand(nrows,1).*DELTA;
        Y=A(:,2*j)./(A(:,2*j)+A(:,2*j-1));
        plot(X,Y,'.');
    endfor
endfor
xlabel('SVLEN'); xticks([1:N_LENGTHS]); xticklabels(LABELS); title('INS inside TR / Mendelian error'); axis square; grid on; set(gca,'fontsize',FONT_SIZE);

%
%
%
% subplot(2,3,2);
% A=load(sprintf('%s/nr_il_mendelian.csv',SOURCE_DIR));
% imagesc(A.'); colormap(COLORMAP); colorbar; axis square;
% xlabel('n repeats in ref'); xticks([1:length(NR_IDS)]); xticklabels(NR_IDS);
% ylabel('interval length'); yticks([1:length(IL_IDS)]); yticklabels(IL_IDS);
% title('Mendelian error (avg over 5 trios)');
%
% subplot(2,3,3);
% A=load(sprintf('%s/ml_nr_mendelian.csv',SOURCE_DIR));
% imagesc(A); colormap(COLORMAP); colorbar; axis square;
% ylabel('motif length'); yticks([1:length(ML_IDS)]); yticklabels(ML_IDS);
% xlabel('n repeats in ref'); xticks([1:length(NR_IDS)]); xticklabels(NR_IDS);
% title('Mendelian error (avg over 5 trios)');
%
%
%
%
%
%
%
%
% figure(2);
%
% subplot(2,3,4);
% A=load(sprintf('%s/ml_il_counts.csv',SOURCE_DIR));
% imagesc((A.')./TOTAL_RECORDS_07_50BP); colormap('hot'); colorbar; axis square;
% xlabel('motif length'); xticks([1:length(ML_IDS)]); xticklabels(ML_IDS);
% ylabel('interval length'); yticks([1:length(IL_IDS)]); yticklabels(IL_IDS);
% title('Fraction of all calls \geq 50bp');
%
% subplot(2,3,5);
% A=load(sprintf('%s/nr_il_counts.csv',SOURCE_DIR));
% imagesc((A.')./TOTAL_RECORDS_07_50BP); colormap('hot'); colorbar; axis square;
% xlabel('n repeats in ref'); xticks([1:length(NR_IDS)]); xticklabels(NR_IDS);
% ylabel('interval length'); yticks([1:length(IL_IDS)]); yticklabels(IL_IDS);
% title('Fraction of all calls \geq 50bp');
%
% subplot(2,3,6);
% A=load(sprintf('%s/ml_nr_counts.csv',SOURCE_DIR));
% imagesc(A./TOTAL_RECORDS_07_50BP); colormap('hot'); colorbar; axis square;
% ylabel('motif length'); yticks([1:length(ML_IDS)]); yticklabels(ML_IDS);
% xlabel('n repeats in ref'); xticks([1:length(NR_IDS)]); xticklabels(NR_IDS);
% title('Fraction of all calls \geq 50bp');




