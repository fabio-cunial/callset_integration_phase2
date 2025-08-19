SOURCE_DIR='.';
FONT_SIZE=18;
DELTA=0.4;

figure(1);




A=load('trios_all.csv');
[nrows,ncolumns]=size(A);
% Columns of A:
% ${KANPIG_N_GOOD_ALT},${KANPIG_N_MERR},
% ${F07_N_GOOD_ALT},${F07_N_MERR},
% ${F09_N_GOOD_ALT},${F09_N_MERR},
% ${COHORT_MERGED_07_N_GOOD_ALT},${COHORT_MERGED_07_N_MERR},
% ${COHORT_MERGED_09_N_GOOD_ALT},${COHORT_MERGED_09_N_MERR},
% ${COHORT_REGENOTYPED_07_N_GOOD_ALT},${COHORT_REGENOTYPED_07_N_MERR},
% ${COHORT_REGENOTYPED_09_N_GOOD_ALT},${COHORT_REGENOTYPED_09_N_MERR}

# All, 07.
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

# All, 09.
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




A=load('trios_tr.csv');
[nrows,ncolumns]=size(A);

# Inside TRs, 07.
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

# Inside TRs, 09.
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




A=load('trios_not_tr.csv');
[nrows,ncolumns]=size(A);

# Outside TRs, 07.
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

# Outside TRs, 09.
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
