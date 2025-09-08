SOURCE_DIR='/Users/fcunial/Downloads/plotDenovoMatrix';
FONT_SIZE=14;
SAMPLE_ID='1004931';
TYPE='all';


% 0/1
figure(1);
subplot(2,2,1);
A=load(sprintf('%s/%s_cohort_regenotyped_07_%s_gtmatrix.txt_denovo_01_adAlt.txt',SOURCE_DIR,SAMPLE_ID,TYPE));
imagesc(log10(A)); colormap('hot'); axis square; ylabel('AD\_CHILD'); xlabel('AD\_FATHER+AD\_MOTHER'); title('De novos in the child with GT=0/1 -- AD ALT'); set(gca,'fontsize',FONT_SIZE);

subplot(2,2,2);
A=load(sprintf('%s/%s_cohort_regenotyped_07_all_gtmatrix.txt_nondenovo_01_adAlt.txt',SOURCE_DIR,SAMPLE_ID));
imagesc(log10(A)); colormap('hot'); axis square; ylabel('AD\_CHILD'); xlabel('AD\_FATHER+AD\_MOTHER'); title('Not de novos occurring in the child with GT=0/1 -- AD ALT'); set(gca,'fontsize',FONT_SIZE);


subplot(2,2,3);
A=load(sprintf('%s/%s_cohort_regenotyped_07_all_gtmatrix.txt_denovo_01_adRef.txt',SOURCE_DIR,SAMPLE_ID));
imagesc(log10(A)); colormap('hot'); axis square; ylabel('AD\_CHILD'); xlabel('AD\_FATHER+AD\_MOTHER'); title('De novos in the child with GT=0/1 -- AD REF'); set(gca,'fontsize',FONT_SIZE);

subplot(2,2,4);
A=load(sprintf('%s/%s_cohort_regenotyped_07_all_gtmatrix.txt_nondenovo_01_adRef.txt',SOURCE_DIR,SAMPLE_ID));
imagesc(log10(A)); colormap('hot'); axis square; ylabel('AD\_CHILD'); xlabel('AD\_FATHER+AD\_MOTHER'); title('Not de novos occurring in the child with GT=0/1 -- AD REF'); set(gca,'fontsize',FONT_SIZE);


% 1/1
figure(2);
subplot(2,2,1);
A=load(sprintf('%s/%s_cohort_regenotyped_07_%s_gtmatrix.txt_denovo_11_adAlt.txt',SOURCE_DIR,SAMPLE_ID,TYPE));
imagesc(log10(A)); colormap('hot'); axis square; ylabel('AD\_CHILD'); xlabel('AD\_FATHER+AD\_MOTHER'); title('De novos in the child with GT=1/1 -- AD ALT'); set(gca,'fontsize',FONT_SIZE);

subplot(2,2,2);
A=load(sprintf('%s/%s_cohort_regenotyped_07_all_gtmatrix.txt_nondenovo_11_adAlt.txt',SOURCE_DIR,SAMPLE_ID));
imagesc(log10(A)); colormap('hot'); axis square; ylabel('AD\_CHILD'); xlabel('AD\_FATHER+AD\_MOTHER'); title('Not de novos occurring in the child with GT=1/1 -- AD ALT'); set(gca,'fontsize',FONT_SIZE);


subplot(2,2,3);
A=load(sprintf('%s/%s_cohort_regenotyped_07_all_gtmatrix.txt_denovo_11_adRef.txt',SOURCE_DIR,SAMPLE_ID));
imagesc(log10(A)); colormap('hot'); axis square; ylabel('AD\_CHILD'); xlabel('AD\_FATHER+AD\_MOTHER'); title('De novos in the child with GT=1/1 -- AD REF'); set(gca,'fontsize',FONT_SIZE);

subplot(2,2,4);
A=load(sprintf('%s/%s_cohort_regenotyped_07_all_gtmatrix.txt_nondenovo_11_adRef.txt',SOURCE_DIR,SAMPLE_ID));
imagesc(log10(A)); colormap('hot'); axis square; ylabel('AD\_CHILD'); xlabel('AD\_FATHER+AD\_MOTHER'); title('Not de novos occurring in the child with GT=1/1 -- AD REF'); set(gca,'fontsize',FONT_SIZE);
















%
%
%
% SOURCE_DIR='.';
% FONT_SIZE=12;
% DELTA=0.4;
% %LABELS={'V1','ND1','ND2', 'ND3','HS1','HS2', 'HS3','HS4','SZ1', 'SZ2','SS1','SS2', 'MQ1','MQ2','MQ3', 'MQ4'};
% %LABELS={'V1','SZ1','SZ2','SS1','SS2','MQ1','MQ2','MQ3','MQ4'};
% LABELS={'V1','FX2'};
% N_ARGS=length(LABELS);
% %N_TOTAL_ARGS=16;
% N_TOTAL_ARGS=2;
% MAX_MENDELIAN_ERROR=0.15;
% MAX_DE_NOVO=0.16;
%
% % Columns:
% % ${N_GOOD_ALT_i},${N_MERR_i},
% % ...,
% %
% % ${DENOVO_i},
% % ...
%
%
% figure(1);
%
% # All, 07.
% A=load(sprintf('%s/squish_trios_all.csv',SOURCE_DIR));
%
% [nrows,ncolumns]=size(A);
%
% subplot(2,3,1); hold on;
% for i=[1:N_ARGS]
%     X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
%     Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
%     plot(X,Y,'.');
% endfor
% title('All variants, 07.'); axis([0,N_ARGS+1,0,MAX_MENDELIAN_ERROR]); grid on;
% xticks([1:N_ARGS]); xticklabels(LABELS);
% ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);
%
% subplot(2,3,4); hold on;
% for i=[1:N_ARGS]
%     X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
%     Y=A(:,N_TOTAL_ARGS*2+i);
%     plot(X,Y,'.');
% endfor
% title('All variants, 07.'); axis([0,N_ARGS+1,0,MAX_DE_NOVO]); grid on;
% xticks([1:N_ARGS]); xticklabels(LABELS);
% ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);
%
%
% # TR, 07.
% A=load(sprintf('%s/squish_trios_tr.csv',SOURCE_DIR));
% [nrows,ncolumns]=size(A);
%
% subplot(2,3,2); hold on;
% for i=[1:N_ARGS]
%     X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
%     Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
%     plot(X,Y,'.');
% endfor
% title('Inside TR, 07.'); axis([0,N_ARGS+1,0,MAX_MENDELIAN_ERROR]); grid on;
% xticks([1:N_ARGS]); xticklabels(LABELS);
% ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);
%
% subplot(2,3,5); hold on;
% for i=[1:N_ARGS]
%     X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
%     Y=A(:,N_TOTAL_ARGS*2+i);
%     plot(X,Y,'.');
% endfor
% title('Inside TR, 07.'); axis([0,N_ARGS+1,0,MAX_DE_NOVO]); grid on;
% xticks([1:N_ARGS]); xticklabels(LABELS);
% ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);
%
%
% # Not TR, 07.
% A=load(sprintf('%s/squish_trios_not_tr.csv',SOURCE_DIR));
% [nrows,ncolumns]=size(A);
%
% subplot(2,3,3); hold on;
% for i=[1:N_ARGS]
%     X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
%     Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
%     plot(X,Y,'.');
% endfor
% title('Outside TR, 07.'); axis([0,N_ARGS+1,0,MAX_MENDELIAN_ERROR]); grid on;
% xticks([1:N_ARGS]); xticklabels(LABELS);
% ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);
%
% subplot(2,3,6); hold on;
% for i=[1:N_ARGS]
%     X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
%     Y=A(:,N_TOTAL_ARGS*2+i);
%     plot(X,Y,'.');
% endfor
% title('Outside TR, 07.'); axis([0,N_ARGS+1,0,MAX_DE_NOVO]); grid on;
% xticks([1:N_ARGS]); xticklabels(LABELS);
% ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);
