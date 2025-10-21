INPUT_DIR_102='/Users/fcunial/Downloads/BenchCohortSamples_PersonalizedCohortVcf/kanpig_1_0_2';
INPUT_DIR_110='/Users/fcunial/Downloads/BenchCohortSamples_PersonalizedCohortVcf/kanpig_1_1_0';
FONT_SIZE=14;
MIN_N_SAMPLES_110=[2, 3, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
MIN_N_SAMPLES_102=[1, 2, 3, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
LABELS={"v1", "1", "2", "3", "4", "8", "16", "32", "64", "128", "256", "512", "1024", "2048"};
DELTA=0.4;


% All calls
A=load(sprintf('%s/personalized_all.csv',INPUT_DIR_110));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # V1, kanpig 1.0.2.
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    subplot(1,3,1); hold on; plot(X,P,'.b'); plot(X,R,'.r'); plot(X,F,'.g'); plot(X,C,'.m');
    # Personalized, kanpig 1.1.0.
    X=[1:length(MIN_N_SAMPLES_110)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES_110)).*DELTA;
    P=[]; R=[]; F=[]; C=[];
    for j=[1:length(MIN_N_SAMPLES_110)]
        P=[P, A(i,j*4+1)];
        R=[R, A(i,j*4+2)];
        F=[F, A(i,j*4+3)];
        C=[C, A(i,j*4+4)];
    endfor
    subplot(1,3,1); hold on; plot(2+X,P,'ob'); plot(2+X,R,'or'); plot(2+X,F,'og'); plot(2+X,C,'om');
endfor
line ([1 2+length(MIN_N_SAMPLES_110)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 2+length(MIN_N_SAMPLES_110)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
line ([1 2+length(MIN_N_SAMPLES_110)], [A(1,3) A(1,3)], "linestyle", "--", "color", "g");
line ([1 2+length(MIN_N_SAMPLES_110)], [A(1,4) A(1,4)], "linestyle", "--", "color", "m");

A=load(sprintf('%s/personalized_all.csv',INPUT_DIR_102));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # Personalized, kanpig 1.0.2.
    X=[1:length(MIN_N_SAMPLES_102)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES_102)).*DELTA;
    P=[]; R=[]; F=[]; C=[];
    for j=[1:length(MIN_N_SAMPLES_102)]
        P=[P, A(i,j*4+1)];
        R=[R, A(i,j*4+2)];
        F=[F, A(i,j*4+3)];
        C=[C, A(i,j*4+4)];
    endfor
    subplot(1,3,1); hold on; plot(1+X,P,'.b'); plot(1+X,R,'.r'); plot(1+X,F,'.g'); plot(1+X,C,'.m');
endfor
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,3) A(1,3)], "linestyle", "--", "color", "g");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,4) A(1,4)], "linestyle", "--", "color", "m");

xticks([1:length(MIN_N_SAMPLES_102)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('All calls, 07'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; legend('precision 1.0.2','recall 1.0.2','F1 1.0.2','GT concordance 1.0.2','precision 1.1.0','recall 1.1.0','F1 1.1.0','GT concordance 1.1.0','location','southoutside'); set(gca,'fontsize',FONT_SIZE);




% Inside TRs
A=load(sprintf('%s/personalized_tr.csv',INPUT_DIR_110));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # V1, kanpig 1.0.2.
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    subplot(1,3,2); hold on; plot(X,P,'.b'); plot(X,R,'.r'); plot(X,F,'.g'); plot(X,C,'.m');
    # Personalized, kanpig 1.1.0.
    X=[1:length(MIN_N_SAMPLES_110)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES_110)).*DELTA;
    P=[]; R=[]; F=[]; C=[];
    for j=[1:length(MIN_N_SAMPLES_110)]
        P=[P, A(i,j*4+1)];
        R=[R, A(i,j*4+2)];
        F=[F, A(i,j*4+3)];
        C=[C, A(i,j*4+4)];
    endfor
    subplot(1,3,2); hold on; plot(2+X,P,'ob'); plot(2+X,R,'or'); plot(2+X,F,'og'); plot(2+X,C,'om');
endfor
line ([1 2+length(MIN_N_SAMPLES_110)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 2+length(MIN_N_SAMPLES_110)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
line ([1 2+length(MIN_N_SAMPLES_110)], [A(1,3) A(1,3)], "linestyle", "--", "color", "g");
line ([1 2+length(MIN_N_SAMPLES_110)], [A(1,4) A(1,4)], "linestyle", "--", "color", "m");

A=load(sprintf('%s/personalized_tr.csv',INPUT_DIR_102));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # Personalized, kanpig 1.0.2.
    X=[1:length(MIN_N_SAMPLES_102)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES_102)).*DELTA;
    P=[]; R=[]; F=[]; C=[];
    for j=[1:length(MIN_N_SAMPLES_102)]
        P=[P, A(i,j*4+1)];
        R=[R, A(i,j*4+2)];
        F=[F, A(i,j*4+3)];
        C=[C, A(i,j*4+4)];
    endfor
    subplot(1,3,2); hold on; plot(1+X,P,'.b'); plot(1+X,R,'.r'); plot(1+X,F,'.g'); plot(1+X,C,'.m');
endfor
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,3) A(1,3)], "linestyle", "--", "color", "g");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,4) A(1,4)], "linestyle", "--", "color", "m");

xticks([1:length(MIN_N_SAMPLES_102)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('Inside TR, 07'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; legend('precision 1.0.2','recall 1.0.2','F1 1.0.2','GT concordance 1.0.2','precision 1.1.0','recall 1.1.0','F1 1.1.0','GT concordance 1.1.0','location','southoutside'); set(gca,'fontsize',FONT_SIZE);




% Outside TRs
A=load(sprintf('%s/personalized_not_tr.csv',INPUT_DIR_110));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # V1, kanpig 1.0.2.
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    subplot(1,3,3); hold on; plot(X,P,'.b'); plot(X,R,'.r'); plot(X,F,'.g'); plot(X,C,'.m');
    # Personalized, kanpig 1.1.0.
    X=[1:length(MIN_N_SAMPLES_110)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES_110)).*DELTA;
    P=[]; R=[]; F=[]; C=[];
    for j=[1:length(MIN_N_SAMPLES_110)]
        P=[P, A(i,j*4+1)];
        R=[R, A(i,j*4+2)];
        F=[F, A(i,j*4+3)];
        C=[C, A(i,j*4+4)];
    endfor
    subplot(1,3,3); hold on; plot(2+X,P,'ob'); plot(2+X,R,'or'); plot(2+X,F,'og'); plot(2+X,C,'om');
endfor
line ([1 2+length(MIN_N_SAMPLES_110)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 2+length(MIN_N_SAMPLES_110)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
line ([1 2+length(MIN_N_SAMPLES_110)], [A(1,3) A(1,3)], "linestyle", "--", "color", "g");
line ([1 2+length(MIN_N_SAMPLES_110)], [A(1,4) A(1,4)], "linestyle", "--", "color", "m");

A=load(sprintf('%s/personalized_not_tr.csv',INPUT_DIR_102));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # Personalized, kanpig 1.0.2.
    X=[1:length(MIN_N_SAMPLES_102)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES_102)).*DELTA;
    P=[]; R=[]; F=[]; C=[];
    for j=[1:length(MIN_N_SAMPLES_102)]
        P=[P, A(i,j*4+1)];
        R=[R, A(i,j*4+2)];
        F=[F, A(i,j*4+3)];
        C=[C, A(i,j*4+4)];
    endfor
    subplot(1,3,3); hold on; plot(1+X,P,'.b'); plot(1+X,R,'.r'); plot(1+X,F,'.g'); plot(1+X,C,'.m');
endfor
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,3) A(1,3)], "linestyle", "--", "color", "g");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,4) A(1,4)], "linestyle", "--", "color", "m");

xticks([1:length(MIN_N_SAMPLES_102)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('Outside TR, 07'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; legend('precision 1.0.2','recall 1.0.2','F1 1.0.2','GT concordance 1.0.2','precision 1.1.0','recall 1.1.0','F1 1.1.0','GT concordance 1.1.0','location','southoutside'); set(gca,'fontsize',FONT_SIZE);























%
% % Outside TRs
% subplot(1,3,3); hold on;
% A=load(sprintf('%s/personalized_not_tr.csv',INPUT_DIR_110));
% [nrows,ncolumns]=size(A);
% for i=[1:nrows]
%     # V1, kanpig 1.0.2.
%     X=1 -DELTA/2 + rand(1,1).*DELTA;
%     P=A(i,1); R=A(i,2);
%     plot(X,P,'.b'); plot(X,R,'.r');
%     # Personalized, kanpig 1.1.0.
%     X=[1:length(MIN_N_SAMPLES)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES)).*DELTA;
%     P=[]; R=[];
%     for j=[1:length(MIN_N_SAMPLES)]
%         P=[P, A(i,(j-1)*2+3)];
%         R=[R, A(i,(j-1)*2+4)];
%     endfor
%     plot(1+X,P,'ob'); plot(1+X,R,'or');
% endfor
% line ([1 1+length(MIN_N_SAMPLES)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
% line ([1 1+length(MIN_N_SAMPLES)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
%
% A=load(sprintf('%s/personalized_not_tr.csv',INPUT_DIR_102));
% [nrows,ncolumns]=size(A);
% for i=[1:nrows]
%     # Personalized, kanpig 1.0.2.
%     X=[1:length(MIN_N_SAMPLES)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES)).*DELTA;
%     P=[]; R=[];
%     for j=[1:length(MIN_N_SAMPLES)]
%         P=[P, A(i,(j-1)*2+3)];
%         R=[R, A(i,(j-1)*2+4)];
%     endfor
%     plot(1+X,P,'.b'); plot(1+X,R,'.r');
% endfor
% line ([1 1+length(MIN_N_SAMPLES)], [A(1,3) A(1,3)], "linestyle", "--", "color", "b");
% line ([1 1+length(MIN_N_SAMPLES)], [A(1,4) A(1,4)], "linestyle", "--", "color", "r");
% xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); title('Outside TR, 07'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; legend('precision 1.0.2','recall 1.0.2','precision 1.1.0','recall 1.1.0','location','southwest'); set(gca,'fontsize',FONT_SIZE);
