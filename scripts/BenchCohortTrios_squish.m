SOURCE_DIR='.';
FONT_SIZE=16;
DELTA=0.4;
%LABELS={'v1','--squish','--ab 0.20','--maxhom 5','--fpenalty 0.2','--fpenalty 0.05','--gpenalty 0.08','--fnmax 1','--fnmax 6'};
LABELS={'V1','SQ','AB','MH','FP1','FP2','GP','FN1','FN2'};

% Columns:
% ${COHORT_REGENOTYPED_07_N_GOOD_ALT},${COHORT_REGENOTYPED_07_N_MERR},
% ${SQUISH_N_GOOD_ALT},${SQUISH_N_MERR},
% ${AB_N_GOOD_ALT},${AB_N_MERR},
% ${MAXHOM_N_GOOD_ALT},${MAXHOM_N_MERR},
% ${FPENALTY1_N_GOOD_ALT},${FPENALTY1_N_MERR},
% ${FPENALTY2_N_GOOD_ALT},${FPENALTY2_N_MERR},
% ${GPENALTY_N_GOOD_ALT},${GPENALTY_N_MERR},
% ${FNMAX1_N_GOOD_ALT},${FNMAX1_N_MERR},
% ${FNMAX2_N_GOOD_ALT},${FNMAX2_N_MERR},
% ${COHORT_REGENOTYPED_07_DENOVO},  // 19-th column (1-based)
% ${SQUISH_DENOVO},
% ${AB_DENOVO},
% ${MAXHOM_DENOVO},
% ${FPENALTY1_DENOVO},
% ${FPENALTY2_DENOVO},
% ${GPENALTY_DENOVO},
% ${FNMAX1_DENOVO},
% ${FNMAX2_DENOVO}


figure(1);

# All, 07.
A=load(sprintf('%s/squish_trios_all.csv',SOURCE_DIR));
[nrows,ncolumns]=size(A);

subplot(2,3,1); hold on;
for i=[1:9]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
    plot(X,Y,'.');
endfor
title('All variants, 07.'); axis([0,10,0,0.15]); axis square; grid on; 
xticks([1:9]); xticklabels(LABELS); xtickangle(45);
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

subplot(2,3,4); hold on;
for i=[19:27]
    X=ones(nrows,1).*(i-18) -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,i);
    plot(X,Y,'.');
endfor
title('All variants, 07.'); axis([0,10,0,0.15]); axis square; grid on; 
xticks([1:9]); xticklabels(LABELS);
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);


# TR, 07.
A=load(sprintf('%s/squish_trios_tr.csv',SOURCE_DIR));
[nrows,ncolumns]=size(A);

subplot(2,3,2); hold on;
for i=[1:9]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
    plot(X,Y,'.');
endfor
title('Inside TR, 07.'); axis([0,10,0,0.15]); axis square; grid on; 
xticks([1:9]); xticklabels(LABELS); xtickangle(45);
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

subplot(2,3,5); hold on;
for i=[19:27]
    X=ones(nrows,1).*(i-18) -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,i);
    plot(X,Y,'.');
endfor
title('Inside TR, 07.'); axis([0,10,0,0.15]); axis square; grid on; 
xticks([1:9]); xticklabels(LABELS);
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);


# Not TR, 07.
A=load(sprintf('%s/squish_trios_not_tr.csv',SOURCE_DIR));
[nrows,ncolumns]=size(A);

subplot(2,3,3); hold on;
for i=[1:9]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
    plot(X,Y,'.');
endfor
title('Outside TR, 07.'); axis([0,10,0,0.15]); axis square; grid on; 
xticks([1:9]); xticklabels(LABELS); xtickangle(45);
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

subplot(2,3,6); hold on;
for i=[19:27]
    X=ones(nrows,1).*(i-18) -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,i);
    plot(X,Y,'.');
endfor
title('Outside TR, 07.'); axis([0,10,0,0.15]); axis square; grid on; 
xticks([1:9]); xticklabels(LABELS);
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);
