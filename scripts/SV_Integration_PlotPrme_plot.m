FONT_SIZE=18;
LABELS={'15x', '30x'};
DELTA=0.4;

EVAL_THRESHOLD='50bp';




% ----------------------------- Precision/recall -------------------------------
figure(1);

% All calls
A=load(sprintf('precision_recall_%s_all.csv',EVAL_THRESHOLD));
[nrows,ncolumns]=size(A);
subplot(1,3,1); hold on;
% 15x
for i=[1:5]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
endfor
AVG_PRECISION_15X=mean(A(1:5,1))
AVG_RECALL_15X=mean(A(1:5,2))
AVG_GT_CONC_15X=mean(A(1:5,4))

% 30x
for i=[6:10]
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
endfor
ALL_AVG_PRECISION_30X=mean(A(6:10,1))
ALL_AVG_RECALL_30X=mean(A(6:10,2))
ALL_AVG_GT_CONC_30X=mean(A(6:10,4))

xticks([1:2]); xticklabels(LABELS); title("All records\nControl samples, whole genome."); grid on; axis([0,3,0,1]); axis square; 
legend('Precision','Recall','GT concordance', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);

% Inside TRs
A=load(sprintf('precision_recall_%s_tr.csv',EVAL_THRESHOLD));
[nrows,ncolumns]=size(A);
subplot(1,3,2); hold on; 
% 15x
for i=[1:5]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
endfor
TR_AVG_PRECISION_15X=mean(A(1:5,1))
TR_AVG_RECALL_15X=mean(A(1:5,2))
TR_AVG_GT_CONC_15X=mean(A(1:5,4))

% 30x
for i=[6:10]
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
endfor
TR_AVG_PRECISION_30X=mean(A(6:10,1))
TR_AVG_RECALL_30X=mean(A(6:10,2))
TR_AVG_GT_CONC_30X=mean(A(6:10,4))

xticks([1:2]); xticklabels(LABELS); title("Inside TRs\nControl samples, whole genome."); grid on; axis([0,3,0,1]); axis square; 
legend('Precision','Recall','GT concordance', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);

% Outside TRs
A=load(sprintf('precision_recall_%s_not_tr.csv',EVAL_THRESHOLD));
[nrows,ncolumns]=size(A);
subplot(1,3,3); hold on; 
% 15x
for i=[1:5]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
endfor
NOT_TR_AVG_PRECISION_15X=mean(A(1:5,1))
NOT_TR_AVG_RECALL_15X=mean(A(1:5,2))
NOT_TR_AVG_GT_CONC_15X=mean(A(1:5,4))

% 30x
for i=[6:10]
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
endfor
NOT_TR_AVG_PRECISION_30X=mean(A(6:10,1))
NOT_TR_AVG_RECALL_30X=mean(A(6:10,2))
NOT_TR_AVG_GT_CONC_30X=mean(A(6:10,4))

xticks([1:2]); xticklabels(LABELS); title("Outside TRs\nControl samples, whole genome."); grid on; axis([0,3,0,1]); axis square; 
legend('Precision','Recall','GT concordance', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);






% ----------------------------- Mendelian error --------------------------------
SUFFIX='_no_missing';

figure(2);
LABELS={"15x\n controls","15x\n AoU","30x\n AoU"};

% All records
B=load(sprintf('mendelian_error_%s_all%s.csv',EVAL_THRESHOLD,SUFFIX));
subplot(1,3,1); hold on;
[nrows,ncolumns]=size(B);
for i=[1:3]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
endfor
ALL_AVG_ME_CONTROLS_15X=mean(B(1:3,2)./(B(1:3,1)+B(1:3,2)))
for i=[4:7]
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
endfor
ALL_AVG_ME_AOU_15X=mean(B(4:7,2)./(B(4:7,1)+B(4:7,2)))
for i=[8:10]
    X=3 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
endfor
ALL_AVG_ME_AOU_30X=mean(B(8:10,2)./(B(8:10,1)+B(8:10,2)))
ylabel('Mendelian error rate'); xticks([1:3]); xticklabels(LABELS); title('All records, whole genome.'); grid on; axis([0,4,0,0.15]); axis square; set(gca,'fontsize',FONT_SIZE);

% Inside TRs
B=load(sprintf('mendelian_error_%s_tr%s.csv',EVAL_THRESHOLD,SUFFIX));
subplot(1,3,2); hold on;
[nrows,ncolumns]=size(B);
for i=[1:3]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
endfor
TR_AVG_ME_CONTROLS_15X=mean(B(1:3,2)./(B(1:3,1)+B(1:3,2)))
for i=[4:6]
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
endfor
TR_AVG_ME_AOU_15X=mean(B(4:6,2)./(B(4:6,1)+B(4:6,2)))
for i=[7:9]
    X=3 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
endfor
TR_AVG_ME_AOU_30X=mean(B(7:9,2)./(B(7:9,1)+B(7:9,2)))
ylabel('Mendelian error rate'); xticks([1:3]); xticklabels(LABELS); title('Inside TRs, whole genome.'); grid on; axis([0,4,0,0.15]); axis square; set(gca,'fontsize',FONT_SIZE);

% Outside TRs
B=load(sprintf('mendelian_error_%s_not_tr%s.csv',EVAL_THRESHOLD,SUFFIX));
subplot(1,3,3); hold on;
[nrows,ncolumns]=size(B);
for i=[1:3]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
endfor
NOT_TR_AVG_ME_CONTROLS_15X=mean(B(1:3,2)./(B(1:3,1)+B(1:3,2)))
for i=[4:6]
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
endfor
NOT_TR_AVG_ME_AOU_15X=mean(B(4:6,2)./(B(4:6,1)+B(4:6,2)))
for i=[7:9]
    X=3 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,2)./(B(i,1)+B(i,2)); plot(X,Y,'.b');
endfor
NOT_TR_AVG_ME_AOU_30X=mean(B(7:9,2)./(B(7:9,1)+B(7:9,2)))
ylabel('Mendelian error rate'); xticks([1:3]); xticklabels(LABELS); title('Outside TRs, whole genome.'); grid on; axis([0,4,0,0.15]); axis square; set(gca,'fontsize',FONT_SIZE);





% ------------------------------ De novo rate ----------------------------------
figure(3);

% All records
B=load(sprintf('denovo_%s_all%s.csv',EVAL_THRESHOLD,SUFFIX));
subplot(1,3,1); hold on;
[nrows,ncolumns]=size(B);
for i=[1:3]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,1)./B(i,2); plot(X,Y,'.b');
endfor
ALL_AVG_DENOVO_CONTROLS_15X=mean(B(1:3,1)./B(1:3,2))
for i=[4:6]
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,1)./B(i,2); plot(X,Y,'.b');
endfor
ALL_AVG_DENOVO_AOU_15X=mean(B(4:6,1)./B(4:6,2))
for i=[7:9]
    X=3 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,1)./B(i,2); plot(X,Y,'.b');
endfor
ALL_AVG_DENOVO_AOU_30X=mean(B(7:9,1)./B(7:9,2))
ylabel('De novo rate'); xticks([1:3]); xticklabels(LABELS); title('All records, whole genome.'); grid on; axis([0,4,0,0.15]); axis square; set(gca,'fontsize',FONT_SIZE);

% Inside TRs
B=load(sprintf('denovo_%s_tr%s.csv',EVAL_THRESHOLD,SUFFIX));
subplot(1,3,2); hold on;
[nrows,ncolumns]=size(B);
for i=[1:3]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,1)./B(i,2); plot(X,Y,'.b');
endfor
TR_AVG_DENOVO_CONTROLS_15X=mean(B(1:3,1)./B(1:3,2))
for i=[4:6]
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,1)./B(i,2); plot(X,Y,'.b');
endfor
TR_AVG_DENOVO_AOU_15X=mean(B(4:6,1)./B(4:6,2))
for i=[7:9]
    X=3 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,1)./B(i,2); plot(X,Y,'.b');
endfor
TR_AVG_DENOVO_AOU_30X=mean(B(7:9,1)./B(7:9,2))
ylabel('De novo rate'); xticks([1:3]); xticklabels(LABELS); title('Inside TRs, whole genome.'); grid on; axis([0,4,0,0.15]); axis square; set(gca,'fontsize',FONT_SIZE);

% Outside TRs
B=load(sprintf('denovo_%s_not_tr%s.csv',EVAL_THRESHOLD,SUFFIX));
subplot(1,3,3); hold on;
[nrows,ncolumns]=size(B);
for i=[1:3]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,1)./B(i,2); plot(X,Y,'.b');
endfor
NOT_TR_AVG_DENOVO_CONTROLS_15X=mean(B(1:3,1)./B(1:3,2))
for i=[4:6]
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,1)./B(i,2); plot(X,Y,'.b');
endfor
NOT_TR_AVG_DENOVO_AOU_15X=mean(B(4:6,1)./B(4:6,2))
for i=[7:9]
    X=3 -DELTA/2 + rand(1,1).*DELTA;
    Y=B(i,1)./B(i,2); plot(X,Y,'.b');
endfor
NOT_TR_AVG_DENOVO_AOU_30X=mean(B(7:9,1)./B(7:9,2))
ylabel('De novo rate'); xticks([1:3]); xticklabels(LABELS); title('Outside TRs, whole genome.'); grid on; axis([0,4,0,0.15]); axis square; set(gca,'fontsize',FONT_SIZE);
