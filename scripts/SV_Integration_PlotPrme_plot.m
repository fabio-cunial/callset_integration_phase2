EVAL_THRESHOLD='20bp_49bp';  # 50bp_10000bp    20bp_49bp

N_SAMPLES_PRECISION_RECALL_15X=40;
N_SAMPLES_PRECISION_RECALL_30X=10;
N_TRIOS_15X_CONTROL=0; #3;
N_TRIOS_15X_AOU=0; #7;
N_TRIOS_30X_AOU=0; #5;

FONT_SIZE=18;
DELTA=0.4;




% ----------------------------- Precision/recall -------------------------------
figure(1);
LABELS={'15x', '30x'};

% All calls
A=load(sprintf('precision_recall_%s_all.csv',EVAL_THRESHOLD));
[nrows,ncolumns]=size(A);
subplot(1,3,1); hold on;
% 15x
for i=[1:N_SAMPLES_PRECISION_RECALL_15X]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
endfor
ALL_AVG_PRECISION_15X=mean(A(1:N_SAMPLES_PRECISION_RECALL_15X,1));
ALL_AVG_RECALL_15X=mean(A(1:N_SAMPLES_PRECISION_RECALL_15X,2));
ALL_AVG_GT_CONC_15X=mean(A(1:N_SAMPLES_PRECISION_RECALL_15X,4));

% 30x
for i=[N_SAMPLES_PRECISION_RECALL_15X+1:N_SAMPLES_PRECISION_RECALL_15X+N_SAMPLES_PRECISION_RECALL_30X]
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
endfor
ALL_AVG_PRECISION_30X=mean(A(N_SAMPLES_PRECISION_RECALL_15X+1:N_SAMPLES_PRECISION_RECALL_15X+N_SAMPLES_PRECISION_RECALL_30X,1));
ALL_AVG_RECALL_30X=mean(A(N_SAMPLES_PRECISION_RECALL_15X+1:N_SAMPLES_PRECISION_RECALL_15X+N_SAMPLES_PRECISION_RECALL_30X,2));
ALL_AVG_GT_CONC_30X=mean(A(N_SAMPLES_PRECISION_RECALL_15X+1:N_SAMPLES_PRECISION_RECALL_15X+N_SAMPLES_PRECISION_RECALL_30X,4));

xticks([1:2]); xticklabels(LABELS); title("All records\nControl samples, whole genome."); grid on; axis([0,3,0,1]); axis square; 
legend('Precision','Recall','GT concordance', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);

% Inside TRs
A=load(sprintf('precision_recall_%s_tr.csv',EVAL_THRESHOLD));
[nrows,ncolumns]=size(A);
subplot(1,3,2); hold on; 
% 15x
for i=[1:N_SAMPLES_PRECISION_RECALL_15X]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
endfor
TR_AVG_PRECISION_15X=mean(A(1:N_SAMPLES_PRECISION_RECALL_15X,1));
TR_AVG_RECALL_15X=mean(A(1:N_SAMPLES_PRECISION_RECALL_15X,2));
TR_AVG_GT_CONC_15X=mean(A(1:N_SAMPLES_PRECISION_RECALL_15X,4));

% 30x
for i=[N_SAMPLES_PRECISION_RECALL_15X+1:N_SAMPLES_PRECISION_RECALL_15X+N_SAMPLES_PRECISION_RECALL_30X]
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
endfor
TR_AVG_PRECISION_30X=mean(A(N_SAMPLES_PRECISION_RECALL_15X+1:N_SAMPLES_PRECISION_RECALL_15X+N_SAMPLES_PRECISION_RECALL_30X,1));
TR_AVG_RECALL_30X=mean(A(N_SAMPLES_PRECISION_RECALL_15X+1:N_SAMPLES_PRECISION_RECALL_15X+N_SAMPLES_PRECISION_RECALL_30X,2));
TR_AVG_GT_CONC_30X=mean(A(N_SAMPLES_PRECISION_RECALL_15X+1:N_SAMPLES_PRECISION_RECALL_15X+N_SAMPLES_PRECISION_RECALL_30X,4));

xticks([1:2]); xticklabels(LABELS); title("Inside TRs\nControl samples, whole genome."); grid on; axis([0,3,0,1]); axis square; 
legend('Precision','Recall','GT concordance', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);

% Outside TRs
A=load(sprintf('precision_recall_%s_not_tr.csv',EVAL_THRESHOLD));
[nrows,ncolumns]=size(A);
subplot(1,3,3); hold on; 
% 15x
for i=[1:N_SAMPLES_PRECISION_RECALL_15X]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
endfor
NOT_TR_AVG_PRECISION_15X=mean(A(1:N_SAMPLES_PRECISION_RECALL_15X,1));
NOT_TR_AVG_RECALL_15X=mean(A(1:N_SAMPLES_PRECISION_RECALL_15X,2));
NOT_TR_AVG_GT_CONC_15X=mean(A(1:N_SAMPLES_PRECISION_RECALL_15X,4));

% 30x
for i=[N_SAMPLES_PRECISION_RECALL_15X+1:N_SAMPLES_PRECISION_RECALL_15X+N_SAMPLES_PRECISION_RECALL_30X]
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
endfor
NOT_TR_AVG_PRECISION_30X=mean(A(N_SAMPLES_PRECISION_RECALL_15X+1:N_SAMPLES_PRECISION_RECALL_15X+N_SAMPLES_PRECISION_RECALL_30X,1));
NOT_TR_AVG_RECALL_30X=mean(A(N_SAMPLES_PRECISION_RECALL_15X+1:N_SAMPLES_PRECISION_RECALL_15X+N_SAMPLES_PRECISION_RECALL_30X,2));
NOT_TR_AVG_GT_CONC_30X=mean(A(N_SAMPLES_PRECISION_RECALL_15X+1:N_SAMPLES_PRECISION_RECALL_15X+N_SAMPLES_PRECISION_RECALL_30X,4));

xticks([1:2]); xticklabels(LABELS); title("Outside TRs\nControl samples, whole genome."); grid on; axis([0,3,0,1]); axis square; 
legend('Precision','Recall','GT concordance', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);






% ----------------------------- Mendelian error --------------------------------
SUFFIX='_no_missing';

if (N_TRIOS_15X_CONTROL > 0)
    figure(2);
    LABELS={"15x\n controls","15x\n AoU","30x\n AoU"};

    % All records
    B=load(sprintf('mendelian_error_%s_all%s.csv',EVAL_THRESHOLD,SUFFIX));
    subplot(1,3,1); hold on;
    [nrows,ncolumns]=size(B);
    for i=[1:N_TRIOS_15X_CONTROL]
        X=1 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
    endfor
    for i=[N_TRIOS_15X_CONTROL+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU]
        X=2 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
    endfor
    ALL_AVG_ME_15X=mean(B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,2)./(B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,1)+B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,2)));
    for i=[N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU]
        X=3 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
    endfor
    ALL_AVG_ME_30X=mean(B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,2)./(B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,1)+B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,2)));
    ylabel('Mendelian error rate'); xticks([1:3]); xticklabels(LABELS); title('All records, whole genome.'); grid on; axis([0,4,0,0.15]); axis square; set(gca,'fontsize',FONT_SIZE);

    % Inside TRs
    B=load(sprintf('mendelian_error_%s_tr%s.csv',EVAL_THRESHOLD,SUFFIX));
    subplot(1,3,2); hold on;
    [nrows,ncolumns]=size(B);
    for i=[1:N_TRIOS_15X_CONTROL]
        X=1 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
    endfor
    for i=[N_TRIOS_15X_CONTROL+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU]
        X=2 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
    endfor
    TR_AVG_ME_15X=mean(B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,2)./(B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,1)+B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,2)));
    for i=[N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU]
        X=3 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
    endfor
    TR_AVG_ME_30X=mean(B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,2)./(B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,1)+B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,2)));
    ylabel('Mendelian error rate'); xticks([1:3]); xticklabels(LABELS); title('Inside TRs, whole genome.'); grid on; axis([0,4,0,0.15]); axis square; set(gca,'fontsize',FONT_SIZE);

    % Outside TRs
    B=load(sprintf('mendelian_error_%s_not_tr%s.csv',EVAL_THRESHOLD,SUFFIX));
    subplot(1,3,3); hold on;
    [nrows,ncolumns]=size(B);
    for i=[1:N_TRIOS_15X_CONTROL]
        X=1 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
    endfor
    for i=[N_TRIOS_15X_CONTROL+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU]
        X=2 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
    endfor
    NOT_TR_AVG_ME_15X=mean(B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,2)./(B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,1)+B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,2)));
    for i=[N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU]
        X=3 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
    endfor
    NOT_TR_AVG_ME_30X=mean(B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,2)./(B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,1)+B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,2)));
    ylabel('Mendelian error rate'); xticks([1:3]); xticklabels(LABELS); title('Outside TRs, whole genome.'); grid on; axis([0,4,0,0.15]); axis square; set(gca,'fontsize',FONT_SIZE);
else
    ALL_AVG_ME_15X=0;
    ALL_AVG_ME_30X=0;
    TR_AVG_ME_15X=0;
    TR_AVG_ME_30X=0;
    NOT_TR_AVG_ME_15X=0;
    NOT_TR_AVG_ME_30X=0;
endif




% ------------------------------ De novo rate ----------------------------------
SUFFIX='_no_missing';

if (N_TRIOS_15X_CONTROL > 0)
    figure(3);

    % All records
    B=load(sprintf('denovo_%s_all%s.csv',EVAL_THRESHOLD,SUFFIX));
    subplot(1,3,1); hold on;
    [nrows,ncolumns]=size(B);
    for i=[1:N_TRIOS_15X_CONTROL]
        X=1 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,1)./B(i,2); plot(X,Y,'.b');
    endfor
    for i=[N_TRIOS_15X_CONTROL+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU]
        X=2 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,1)./B(i,2); plot(X,Y,'.b');
    endfor
    ALL_AVG_DENOVO_15X=mean(B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,1)./B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,2));
    for i=[N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU]
        X=3 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,1)./B(i,2); plot(X,Y,'.b');
    endfor
    ALL_AVG_DENOVO_30X=mean(B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,1)./B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,2));
    ylabel('De novo rate'); xticks([1:3]); xticklabels(LABELS); title('All records, whole genome.'); grid on; axis([0,4,0,0.15]); axis square; set(gca,'fontsize',FONT_SIZE);

    % Inside TRs
    B=load(sprintf('denovo_%s_tr%s.csv',EVAL_THRESHOLD,SUFFIX));
    subplot(1,3,2); hold on;
    [nrows,ncolumns]=size(B);
    for i=[1:N_TRIOS_15X_CONTROL]
        X=1 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,1)./B(i,2); plot(X,Y,'.b');
    endfor
    for i=[N_TRIOS_15X_CONTROL+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU]
        X=2 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,1)./B(i,2); plot(X,Y,'.b');
    endfor
    TR_AVG_DENOVO_15X=mean(B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,1)./B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,2));
    for i=[N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU]
        X=3 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,1)./B(i,2); plot(X,Y,'.b');
    endfor
    TR_AVG_DENOVO_30X=mean(B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,1)./B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,2));
    ylabel('De novo rate'); xticks([1:3]); xticklabels(LABELS); title('Inside TRs, whole genome.'); grid on; axis([0,4,0,0.15]); axis square; set(gca,'fontsize',FONT_SIZE);

    % Outside TRs
    B=load(sprintf('denovo_%s_not_tr%s.csv',EVAL_THRESHOLD,SUFFIX));
    subplot(1,3,3); hold on;
    [nrows,ncolumns]=size(B);
    for i=[1:N_TRIOS_15X_CONTROL]
        X=1 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,1)./B(i,2); plot(X,Y,'.b');
    endfor
    for i=[N_TRIOS_15X_CONTROL+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU]
        X=2 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,1)./B(i,2); plot(X,Y,'.b');
    endfor
    NOT_TR_AVG_DENOVO_15X=mean(B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,1)./B(1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU,2));
    for i=[N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU]
        X=3 -DELTA/2 + rand(1,1).*DELTA;
        Y=B(i,1)./B(i,2); plot(X,Y,'.b');
    endfor
    NOT_TR_AVG_DENOVO_30X=mean(B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,1)./B(N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+1:N_TRIOS_15X_CONTROL+N_TRIOS_15X_AOU+N_TRIOS_30X_AOU,2));
    ylabel('De novo rate'); xticks([1:3]); xticklabels(LABELS); title('Outside TRs, whole genome.'); grid on; axis([0,4,0,0.15]); axis square; set(gca,'fontsize',FONT_SIZE);
else
    ALL_AVG_DENOVO_15X=0;
    ALL_AVG_DENOVO_30X=0;
    TR_AVG_DENOVO_15X=0;
    TR_AVG_DENOVO_30X=0;
    NOT_TR_AVG_DENOVO_15X=0;
    NOT_TR_AVG_DENOVO_30X=0;
endif




% ------------------------------ Final output ----------------------------------
fprintf('    |        Precision      |         Recall        |    GT concordance     |  Mendelian error rate |      De novo rate     |\n');
fprintf('    |  All  |   TR  | Non-TR|  All  |   TR  | Non-TR|  All  |  TR   | Non-TR|  All  |  TR   | Non-TR|  All  |  TR   | Non-TR|\n');
fprintf('15x | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |\n', ALL_AVG_PRECISION_15X,TR_AVG_PRECISION_15X,NOT_TR_AVG_PRECISION_15X,  ALL_AVG_RECALL_15X,TR_AVG_RECALL_15X,NOT_TR_AVG_RECALL_15X,  ALL_AVG_GT_CONC_15X,TR_AVG_GT_CONC_15X,NOT_TR_AVG_GT_CONC_15X,  ALL_AVG_ME_15X,TR_AVG_ME_15X,NOT_TR_AVG_ME_15X,  ALL_AVG_DENOVO_15X,TR_AVG_DENOVO_15X,NOT_TR_AVG_DENOVO_15X);
fprintf('30x | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |\n', ALL_AVG_PRECISION_30X,TR_AVG_PRECISION_30X,NOT_TR_AVG_PRECISION_30X,  ALL_AVG_RECALL_30X,TR_AVG_RECALL_30X,NOT_TR_AVG_RECALL_30X,  ALL_AVG_GT_CONC_30X,TR_AVG_GT_CONC_30X,NOT_TR_AVG_GT_CONC_30X,  ALL_AVG_ME_30X,TR_AVG_ME_30X,NOT_TR_AVG_ME_30X,  ALL_AVG_DENOVO_30X,TR_AVG_DENOVO_30X,NOT_TR_AVG_DENOVO_30X);
