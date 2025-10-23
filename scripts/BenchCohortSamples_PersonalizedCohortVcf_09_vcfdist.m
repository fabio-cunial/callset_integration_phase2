INPUT_DIR_102='/Users/fcunial/Downloads/BenchCohortSamples_PersonalizedCohortVcf/09_kanpig_1_0_2_vcfdist';
FONT_SIZE=14;
MIN_N_SAMPLES_102=[2, 4, 8, 16, 32, 64, 128, 256, 512];
LABELS={"v1", "2", "4", "8", "16", "32", "64", "128", "256", "512"};
DELTA=0.4;


% All calls
A=load(sprintf('%s/personalized_all.csv',INPUT_DIR_102));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # V1, kanpig 1.0.2.
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3);
    subplot(1,3,1); hold on; plot(X,P,'.b'); plot(X,R,'.r'); plot(X,F,'.g');
    # Personalized, kanpig 1.0.2.
    X=[1:length(MIN_N_SAMPLES_102)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES_102)).*DELTA;
    P=[]; R=[]; F=[];
    for j=[1:length(MIN_N_SAMPLES_102)]
        P=[P, A(i,j*3+1)];
        R=[R, A(i,j*3+2)];
        F=[F, A(i,j*3+3)];
    endfor
    subplot(1,3,1); hold on; plot(1+X,P,'.b'); plot(1+X,R,'.r'); plot(1+X,F,'.g');
endfor
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,3) A(1,3)], "linestyle", "--", "color", "g");
xticks([1:length(MIN_N_SAMPLES_102)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('All calls, 09'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; legend('precision 1.0.2','recall 1.0.2','F1 1.0.2','location','southoutside'); set(gca,'fontsize',FONT_SIZE);




% Inside TRs
A=load(sprintf('%s/personalized_tr.csv',INPUT_DIR_102));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # V1, kanpig 1.0.2.
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3);
    subplot(1,3,2); hold on; plot(X,P,'.b'); plot(X,R,'.r'); plot(X,F,'.g');
    # Personalized, kanpig 1.0.2.
    X=[1:length(MIN_N_SAMPLES_102)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES_102)).*DELTA;
    P=[]; R=[]; F=[];
    for j=[1:length(MIN_N_SAMPLES_102)]
        P=[P, A(i,j*3+1)];
        R=[R, A(i,j*3+2)];
        F=[F, A(i,j*3+3)];
    endfor
    subplot(1,3,2); hold on; plot(1+X,P,'.b'); plot(1+X,R,'.r'); plot(1+X,F,'.g');
endfor
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,3) A(1,3)], "linestyle", "--", "color", "g");
xticks([1:length(MIN_N_SAMPLES_102)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('Inside TR, 09'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; legend('precision 1.0.2','recall 1.0.2','F1 1.0.2','location','southoutside'); set(gca,'fontsize',FONT_SIZE);




% Outside TRs
A=load(sprintf('%s/personalized_not_tr.csv',INPUT_DIR_102));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # V1, kanpig 1.0.2.
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3);
    subplot(1,3,3); hold on; plot(X,P,'.b'); plot(X,R,'.r'); plot(X,F,'.g');
    # Personalized, kanpig 1.0.2.
    X=[1:length(MIN_N_SAMPLES_102)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES_102)).*DELTA;
    P=[]; R=[]; F=[];
    for j=[1:length(MIN_N_SAMPLES_102)]
        P=[P, A(i,j*3+1)];
        R=[R, A(i,j*3+2)];
        F=[F, A(i,j*3+3)];
    endfor
    subplot(1,3,3); hold on; plot(1+X,P,'.b'); plot(1+X,R,'.r'); plot(1+X,F,'.g');
endfor
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
line ([1 1+length(MIN_N_SAMPLES_102)], [A(1,3) A(1,3)], "linestyle", "--", "color", "g");
xticks([1:length(MIN_N_SAMPLES_102)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('Outside TR, 09'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; legend('precision 1.0.2','recall 1.0.2','F1 1.0.2','location','southoutside'); set(gca,'fontsize',FONT_SIZE);
