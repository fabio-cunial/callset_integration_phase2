SOURCE_DIR='.';
FONT_SIZE=18;
MIN_N_SAMPLES=[2, 3, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
LABELS={"v1", "2", "3", "4", "8", "16", "32", "64", "128", "256", "512", "1024", "2048"};
DELTA=0.4;


% All calls
subplot(1,3,1); hold on;
A=load('personalized_all.csv');
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    X=[1:1+length(MIN_N_SAMPLES)] -DELTA/2 + rand(1,1+length(MIN_N_SAMPLES)).*DELTA;
    P=[]; R=[];
    for j=[1:length(MIN_N_SAMPLES)+1]
        P=[P, A(i,(j-1)*2+1)];
        R=[R, A(i,(j-1)*2+2)];
    endfor
    plot(X,P,'.b'); plot(X,R,'.r');
endfor
line ([1 1+length(MIN_N_SAMPLES)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,3) A(1,3)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,4) A(1,4)], "linestyle", "--", "color", "r");
xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); title('All calls, 07'); grid on; axis([1,length(LABELS),0,1]); axis square; legend('precision','recall','location','southwest'); set(gca,'fontsize',FONT_SIZE);

% Inside TRs
subplot(1,3,2); hold on;
A=load('personalized_tr.csv');
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    X=[1:1+length(MIN_N_SAMPLES)] -DELTA/2 + rand(1,1+length(MIN_N_SAMPLES)).*DELTA;
    P=[]; R=[];
    for j=[1:length(MIN_N_SAMPLES)+1]
        P=[P, A(i,(j-1)*2+1)];
        R=[R, A(i,(j-1)*2+2)];
    endfor
    plot(X,P,'.b'); plot(X,R,'.r');
endfor
line ([1 1+length(MIN_N_SAMPLES)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,3) A(1,3)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,4) A(1,4)], "linestyle", "--", "color", "r");
xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); title('Inside TR, 07'); grid on; axis([1,length(LABELS),0,1]); axis square; legend('precision','recall','location','southwest'); set(gca,'fontsize',FONT_SIZE);

% Outside TRs
subplot(1,3,3); hold on;
A=load('personalized_not_tr.csv');
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    X=[1:1+length(MIN_N_SAMPLES)] -DELTA/2 + rand(1,1+length(MIN_N_SAMPLES)).*DELTA;
    P=[]; R=[];
    for j=[1:length(MIN_N_SAMPLES)+1]
        P=[P, A(i,(j-1)*2+1)];
        R=[R, A(i,(j-1)*2+2)];
    endfor
    plot(X,P,'.b'); plot(X,R,'.r');
endfor
line ([1 1+length(MIN_N_SAMPLES)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,3) A(1,3)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,4) A(1,4)], "linestyle", "--", "color", "r");
xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); title('Outside TR, 07'); grid on; axis([1,length(LABELS),0,1]); axis square; legend('precision','recall','location','southwest'); set(gca,'fontsize',FONT_SIZE);




% betabinomial=load('beta_all.csv');
% subplot(1,3,1); hold on;
% P07=[kanpig(:,1),betabinomial(:,1)];
% R07=[kanpig(:,2),betabinomial(:,2)];
% plot(P07(:,1),R07(:,1),'.r');
% plot(P07(:,2),R07(:,2),'.g');
% title('All calls'); axis([0.8,1,0.6,1]); axis square; grid on; xlabel('Precision'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);
%
% # Inside TRs
% kanpig=load('kanpig_tr.csv');
% betabinomial=load('beta_tr.csv');
% subplot(1,3,2); hold on;
% P07=[kanpig(:,1),betabinomial(:,1)];
% R07=[kanpig(:,2),betabinomial(:,2)];
% plot(P07(:,1),R07(:,1),'.r');
% plot(P07(:,2),R07(:,2),'.g');
% title('Inside TR'); axis([0.8,1,0.6,1]); axis square; grid on; xlabel('Precision'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);
%
% # Outside TRs
% kanpig=load('kanpig_not_tr.csv');
% betabinomial=load('beta_not_tr.csv');
% subplot(1,3,3); hold on;
% P07=[kanpig(:,1),betabinomial(:,1)];
% R07=[kanpig(:,2),betabinomial(:,2)];
% plot(P07(:,1),R07(:,1),'.r');
% plot(P07(:,2),R07(:,2),'.g');
% title('Outside TR'); axis([0.8,1,0.6,1]); axis square; grid on; xlabel('Precision'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);
