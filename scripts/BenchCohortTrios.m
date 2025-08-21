SOURCE_DIR='/Users/fcunial/Downloads/BenchCohortTrios_collapse_prime';
FONT_SIZE=18;
DELTA=0.4;

% Columns:
% ${KANPIG_N_GOOD_ALT},${KANPIG_N_MERR},
% ${F07_N_GOOD_ALT},${F07_N_MERR},
% ${F09_N_GOOD_ALT},${F09_N_MERR},
% ${COHORT_MERGED_07_N_GOOD_ALT},${COHORT_MERGED_07_N_MERR},
% ${COHORT_MERGED_09_N_GOOD_ALT},${COHORT_MERGED_09_N_MERR},
% ${COHORT_REGENOTYPED_07_N_GOOD_ALT},${COHORT_REGENOTYPED_07_N_MERR},
% ${COHORT_REGENOTYPED_09_N_GOOD_ALT},${COHORT_REGENOTYPED_09_N_MERR},
% ${KANPIG_DENOVO},   // 15th column (1-based)
% ${F07_DENOVO},
% ${F09_DENOVO},
% ${COHORT_MERGED_07_DENOVO},
% ${COHORT_MERGED_09_DENOVO},
% ${COHORT_REGENOTYPED_07_DENOVO},
% ${COHORT_REGENOTYPED_09_DENOVO}




A=load(sprintf('%s/trios_all.csv',SOURCE_DIR));
[nrows,ncolumns]=size(A);

# All, 07.
figure(1);
subplot(3,2,1); hold on;
X=ones(nrows,1).*1 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,2)./(A(:,2)+A(:,1));
plot(X,Y,'.');

X=ones(nrows,1).*2 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,4)./(A(:,4)+A(:,3));
plot(X,Y,'.');

X=ones(nrows,1).*3 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,8)./(A(:,8)+A(:,7));
plot(X,Y,'.');

X=ones(nrows,1).*4 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,12)./(A(:,12)+A(:,11));
plot(X,Y,'.');

title('All variants, 07.'); axis([0,5,0,0.14]); axis square; grid on; 
xticks([1:4]); xticklabels({sprintf('kanpig\nintra'),'xgb',sprintf('truvari\ninter'),sprintf('kanpig\ninter')});
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

figure(2);
subplot(3,2,1); hold on;
X=ones(nrows,1).*1 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,15);
plot(X,Y,'.');

X=ones(nrows,1).*2 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,16);
plot(X,Y,'.');

X=ones(nrows,1).*3 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,18);
plot(X,Y,'.');

X=ones(nrows,1).*4 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,20);
plot(X,Y,'.');

title('All variants, 07.'); axis([0,5,0,0.14]); axis square; grid on; 
xticks([1:4]); xticklabels({sprintf('kanpig\nintra'),'xgb',sprintf('truvari\ninter'),sprintf('kanpig\ninter')});
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);



# All, 09.
figure(1);
subplot(3,2,2); hold on;
X=ones(nrows,1).*1 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,2)./(A(:,2)+A(:,1));
plot(X,Y,'.');

X=ones(nrows,1).*2 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,6)./(A(:,6)+A(:,5));
plot(X,Y,'.');

X=ones(nrows,1).*3 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,10)./(A(:,10)+A(:,9));
plot(X,Y,'.');

X=ones(nrows,1).*4 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,14)./(A(:,14)+A(:,13));
plot(X,Y,'.');

title('All variants, 09.'); axis([0,5,0,0.14]); axis square; grid on; 
xticks([1:4]); xticklabels({sprintf('kanpig\nintra'),'xgb',sprintf('truvari\ninter'),sprintf('kanpig\ninter')});
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

figure(2);
subplot(3,2,2); hold on;
X=ones(nrows,1).*1 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,15);
plot(X,Y,'.');

X=ones(nrows,1).*2 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,17);
plot(X,Y,'.');

X=ones(nrows,1).*3 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,19);
plot(X,Y,'.');

X=ones(nrows,1).*4 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,21);
plot(X,Y,'.');

title('All variants, 09.'); axis([0,5,0,0.14]); axis square; grid on; 
xticks([1:4]); xticklabels({sprintf('kanpig\nintra'),'xgb',sprintf('truvari\ninter'),sprintf('kanpig\ninter')});
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);




A=load(sprintf('%s/trios_tr.csv',SOURCE_DIR));
[nrows,ncolumns]=size(A);

# Inside TRs, 07.
figure(1);
subplot(3,2,3); hold on;
X=ones(nrows,1).*1 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,2)./(A(:,2)+A(:,1));
plot(X,Y,'.');

X=ones(nrows,1).*2 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,4)./(A(:,4)+A(:,3));
plot(X,Y,'.');

X=ones(nrows,1).*3 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,8)./(A(:,8)+A(:,7));
plot(X,Y,'.');

X=ones(nrows,1).*4 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,12)./(A(:,12)+A(:,11));
plot(X,Y,'.');

title('Inside TRs, 07.'); axis([0,5,0,0.14]); axis square; grid on; 
xticks([1:4]); xticklabels({sprintf('kanpig\nintra'),'xgb',sprintf('truvari\ninter'),sprintf('kanpig\ninter')});
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

figure(2);
subplot(3,2,3); hold on;
X=ones(nrows,1).*1 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,15);
plot(X,Y,'.');

X=ones(nrows,1).*2 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,16);
plot(X,Y,'.');

X=ones(nrows,1).*3 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,18);
plot(X,Y,'.');

X=ones(nrows,1).*4 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,20);
plot(X,Y,'.');

title('Inside TRs, 07.'); axis([0,5,0,0.14]); axis square; grid on; 
xticks([1:4]); xticklabels({sprintf('kanpig\nintra'),'xgb',sprintf('truvari\ninter'),sprintf('kanpig\ninter')});
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);




# Inside TRs, 09.
figure(1);
subplot(3,2,4); hold on;
X=ones(nrows,1).*1 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,2)./(A(:,2)+A(:,1));
plot(X,Y,'.');

X=ones(nrows,1).*2 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,6)./(A(:,6)+A(:,5));
plot(X,Y,'.');

X=ones(nrows,1).*3 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,10)./(A(:,10)+A(:,9));
plot(X,Y,'.');

X=ones(nrows,1).*4 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,14)./(A(:,14)+A(:,13));
plot(X,Y,'.');

title('Inside TRs, 09.'); axis([0,5,0,0.14]); axis square; grid on; 
xticks([1:4]); xticklabels({sprintf('kanpig\nintra'),'xgb',sprintf('truvari\ninter'),sprintf('kanpig\ninter')});
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

figure(2);
subplot(3,2,4); hold on;
X=ones(nrows,1).*1 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,15);
plot(X,Y,'.');

X=ones(nrows,1).*2 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,17);
plot(X,Y,'.');

X=ones(nrows,1).*3 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,19);
plot(X,Y,'.');

X=ones(nrows,1).*4 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,21);
plot(X,Y,'.');

title('Inside TRs, 09.'); axis([0,5,0,0.14]); axis square; grid on; 
xticks([1:4]); xticklabels({sprintf('kanpig\nintra'),'xgb',sprintf('truvari\ninter'),sprintf('kanpig\ninter')});
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);




A=load(sprintf('%s/trios_not_tr.csv',SOURCE_DIR));
[nrows,ncolumns]=size(A);

# Outside TRs, 07.
figure(1);
subplot(3,2,5); hold on;
X=ones(nrows,1).*1 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,2)./(A(:,2)+A(:,1));
plot(X,Y,'.');

X=ones(nrows,1).*2 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,4)./(A(:,4)+A(:,3));
plot(X,Y,'.');

X=ones(nrows,1).*3 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,8)./(A(:,8)+A(:,7));
plot(X,Y,'.');

X=ones(nrows,1).*4 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,12)./(A(:,12)+A(:,11));
plot(X,Y,'.');

title('Outside TRs, 07.'); axis([0,5,0,0.14]); axis square; grid on; 
xticks([1:4]); xticklabels({sprintf('kanpig\nintra'),'xgb',sprintf('truvari\ninter'),sprintf('kanpig\ninter')});
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);


figure(2);
subplot(3,2,5); hold on;
X=ones(nrows,1).*1 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,15);
plot(X,Y,'.');

X=ones(nrows,1).*2 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,16);
plot(X,Y,'.');

X=ones(nrows,1).*3 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,18);
plot(X,Y,'.');

X=ones(nrows,1).*4 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,20);
plot(X,Y,'.');

title('Outside TRs, 07.'); axis([0,5,0,0.14]); axis square; grid on; 
xticks([1:4]); xticklabels({sprintf('kanpig\nintra'),'xgb',sprintf('truvari\ninter'),sprintf('kanpig\ninter')});
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);




# Outside TRs, 09.
figure(1);
subplot(3,2,6); hold on;
X=ones(nrows,1).*1 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,2)./(A(:,2)+A(:,1));
plot(X,Y,'.');

X=ones(nrows,1).*2 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,6)./(A(:,6)+A(:,5));
plot(X,Y,'.');

X=ones(nrows,1).*3 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,10)./(A(:,10)+A(:,9));
plot(X,Y,'.');

X=ones(nrows,1).*4 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,14)./(A(:,14)+A(:,13));
plot(X,Y,'.');

title('Outside TRs, 09.'); axis([0,5,0,0.14]); axis square; grid on; 
xticks([1:4]); xticklabels({sprintf('kanpig\nintra'),'xgb',sprintf('truvari\ninter'),sprintf('kanpig\ninter')});
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);


figure(2);
subplot(3,2,6); hold on;
X=ones(nrows,1).*1 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,15);
plot(X,Y,'.');

X=ones(nrows,1).*2 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,16);
plot(X,Y,'.');

X=ones(nrows,1).*3 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,18);
plot(X,Y,'.');

X=ones(nrows,1).*4 -DELTA/2 + rand(nrows,1).*DELTA;
Y=A(:,20);
plot(X,Y,'.');

title('Outside TRs, 09.'); axis([0,5,0,0.14]); axis square; grid on; 
xticks([1:4]); xticklabels({sprintf('kanpig\nintra'),'xgb',sprintf('truvari\ninter'),sprintf('kanpig\ninter')});
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);
