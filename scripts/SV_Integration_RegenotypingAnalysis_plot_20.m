FONT_SIZE=14;
MIN_N_SAMPLES=[1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
LABELS={"T", "1", "2", "4", "8", "16", "32", "64", "128", "256", "512", "1024", "2048"};
DELTA=0.4;

EVAL_THRESHOLD='20bp';




% ----------------------------- Precision/recall -------------------------------
figure(1);

% All calls
A=load(sprintf('precision_recall_15x_%s_all.csv',EVAL_THRESHOLD));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # Truvari collapse
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    subplot(1,3,1); hold on; plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
    # Re-genotyping, 15x control.
    X=[1:length(MIN_N_SAMPLES)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES)).*DELTA;
    P=[]; R=[]; F=[]; C=[];
    for j=[1:length(MIN_N_SAMPLES)]
        P=[P, A(i,j*4+1)];
        R=[R, A(i,j*4+2)];
        F=[F, A(i,j*4+3)];
        C=[C, A(i,j*4+4)];
    endfor
    subplot(1,3,1); hold on; plot(1+X,P,'.b'); plot(1+X,R,'.r');
	%plot(1+X,F,'.g');
	plot(1+X,C,'.m');
endfor
line ([1 1+length(MIN_N_SAMPLES)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
%line ([1 1+length(MIN_N_SAMPLES)], [A(1,3) A(1,3)], "linestyle", "--", "color", "g");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,4) A(1,4)], "linestyle", "--", "color", "m");

A=load(sprintf('precision_recall_30x_%s_all.csv',EVAL_THRESHOLD));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # Truvari collapse
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    subplot(1,3,1); hold on; plot(X,P,'ob'); plot(X,R,'or'); 
	%plot(X,F,'.g'); 
	plot(X,C,'om');
    # Re-genotyping, 30x control.
    X=[1:length(MIN_N_SAMPLES)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES)).*DELTA;
    P=[]; R=[]; F=[]; C=[];
    for j=[1:length(MIN_N_SAMPLES)]
        P=[P, A(i,j*4+1)];
        R=[R, A(i,j*4+2)];
        F=[F, A(i,j*4+3)];
        C=[C, A(i,j*4+4)];
    endfor
    subplot(1,3,1); hold on; plot(1+X,P,'ob'); plot(1+X,R,'or'); 
	%plot(1+X,F,'og'); 
	plot(1+X,C,'om');
endfor
xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('All records'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; 
legend('15x precision','15x recall','15x GT concordance', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);


% Inside TRs
A=load(sprintf('precision_recall_15x_%s_tr.csv',EVAL_THRESHOLD));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # Truvari collapse
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    subplot(1,3,2); hold on; plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
    # Re-genotyping, 15x control.
    X=[1:length(MIN_N_SAMPLES)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES)).*DELTA;
    P=[]; R=[]; F=[]; C=[];
    for j=[1:length(MIN_N_SAMPLES)]
        P=[P, A(i,j*4+1)];
        R=[R, A(i,j*4+2)];
        F=[F, A(i,j*4+3)];
        C=[C, A(i,j*4+4)];
    endfor
    subplot(1,3,2); hold on; plot(1+X,P,'.b'); plot(1+X,R,'.r'); 
	%plot(1+X,F,'.g'); 
	plot(1+X,C,'.m');
endfor
line ([1 1+length(MIN_N_SAMPLES)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
%line ([1 1+length(MIN_N_SAMPLES)], [A(1,3) A(1,3)], "linestyle", "--", "color", "g");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,4) A(1,4)], "linestyle", "--", "color", "m");

A=load(sprintf('precision_recall_30x_%s_tr.csv',EVAL_THRESHOLD));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # Truvari collapse
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    subplot(1,3,2); hold on; plot(X,P,'ob'); plot(X,R,'or'); 
	%plot(X,F,'.g'); 
	plot(X,C,'om');
    # Re-genotyping, 30x control.
    X=[1:length(MIN_N_SAMPLES)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES)).*DELTA;
    P=[]; R=[]; F=[]; C=[];
    for j=[1:length(MIN_N_SAMPLES)]
        P=[P, A(i,j*4+1)];
        R=[R, A(i,j*4+2)];
        F=[F, A(i,j*4+3)];
        C=[C, A(i,j*4+4)];
    endfor
    subplot(1,3,2); hold on; plot(1+X,P,'ob'); plot(1+X,R,'or'); 
	%plot(1+X,F,'og'); 
	plot(1+X,C,'om');
endfor
xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('Inside TRs'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; 
legend('15x precision','15x recall','15x GT concordance', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);


% Outside TRs
A=load(sprintf('precision_recall_15x_%s_not_tr.csv',EVAL_THRESHOLD));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # Truvari collapse
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    subplot(1,3,3); hold on; plot(X,P,'.b'); plot(X,R,'.r'); 
	%plot(X,F,'.g'); 
	plot(X,C,'.m');
    # Re-genotyping, 15x control.
    X=[1:length(MIN_N_SAMPLES)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES)).*DELTA;
    P=[]; R=[]; F=[]; C=[];
    for j=[1:length(MIN_N_SAMPLES)]
        P=[P, A(i,j*4+1)];
        R=[R, A(i,j*4+2)];
        F=[F, A(i,j*4+3)];
        C=[C, A(i,j*4+4)];
    endfor
    subplot(1,3,3); hold on; plot(1+X,P,'.b'); plot(1+X,R,'.r'); 
	%plot(1+X,F,'.g'); 
	plot(1+X,C,'.m');
endfor
line ([1 1+length(MIN_N_SAMPLES)], [A(1,1) A(1,1)], "linestyle", "--", "color", "b");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,2) A(1,2)], "linestyle", "--", "color", "r");
%line ([1 1+length(MIN_N_SAMPLES)], [A(1,3) A(1,3)], "linestyle", "--", "color", "g");
line ([1 1+length(MIN_N_SAMPLES)], [A(1,4) A(1,4)], "linestyle", "--", "color", "m");

A=load(sprintf('precision_recall_30x_%s_not_tr.csv',EVAL_THRESHOLD));
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    # Truvari collapse
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2); F=A(i,3); C=A(i,4);
    subplot(1,3,3); hold on; plot(X,P,'ob'); plot(X,R,'or'); 
	%plot(X,F,'.g'); 
	plot(X,C,'om');
    # Re-genotyping, 30x control.
    X=[1:length(MIN_N_SAMPLES)] -DELTA/2 + rand(1,length(MIN_N_SAMPLES)).*DELTA;
    P=[]; R=[]; F=[]; C=[];
    for j=[1:length(MIN_N_SAMPLES)]
        P=[P, A(i,j*4+1)];
        R=[R, A(i,j*4+2)];
        F=[F, A(i,j*4+3)];
        C=[C, A(i,j*4+4)];
    endfor
    subplot(1,3,3); hold on; plot(1+X,P,'ob'); plot(1+X,R,'or'); 
	%plot(1+X,F,'og'); 
	plot(1+X,C,'om');
endfor
xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('Outside TRs'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; 
legend('15x precision','15x recall','15x GT concordance', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);




% ----------------------------- Mendelian error --------------------------------
figure(2);

% All records
A=load(sprintf('mendelian_error_15x_%s_control_all.csv',EVAL_THRESHOLD));
B=load(sprintf('mendelian_error_15x_%s_aou_all.csv',EVAL_THRESHOLD));
C=load(sprintf('mendelian_error_30x_%s_aou_all.csv',EVAL_THRESHOLD));
subplot(1,3,1); hold on;
[nrows,ncolumns]=size(A);
for i=[2:2:ncolumns]
	[nrows,ncolumns]=size(A); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=A(:,i)./(A(:,i)+A(:,i-1)); plot(X,Y,'.r');
	[nrows,ncolumns]=size(B); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=B(:,i)./(B(:,i)+B(:,i-1)); plot(X,Y,'.b');
	[nrows,ncolumns]=size(C); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=C(:,i)./(C(:,i)+C(:,i-1)); plot(X,Y,'ob');
endfor
ylabel('Mendelian error rate'); xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('All records'); grid on; axis([0,length(LABELS)+1,0,0.25]); axis square; set(gca,'fontsize',FONT_SIZE);
legend('15x control','15x AoU','30x AoU', 'location','southoutside'); 


% Inside TRs
A=load(sprintf('mendelian_error_15x_%s_control_tr.csv',EVAL_THRESHOLD));
B=load(sprintf('mendelian_error_15x_%s_aou_tr.csv',EVAL_THRESHOLD));
C=load(sprintf('mendelian_error_30x_%s_aou_tr.csv',EVAL_THRESHOLD));
subplot(1,3,2); hold on;
[nrows,ncolumns]=size(A);
for i=[2:2:ncolumns]
	[nrows,ncolumns]=size(A); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=A(:,i)./(A(:,i)+A(:,i-1)); plot(X,Y,'.r');
	[nrows,ncolumns]=size(B); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=B(:,i)./(B(:,i)+B(:,i-1)); plot(X,Y,'.b');
	[nrows,ncolumns]=size(C); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=C(:,i)./(C(:,i)+C(:,i-1)); plot(X,Y,'ob');
endfor
ylabel('Mendelian error rate'); xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('Inside TRs'); grid on; axis([0,length(LABELS)+1,0,0.25]); axis square; set(gca,'fontsize',FONT_SIZE);
legend('15x control','15x AoU','30x AoU', 'location','southoutside'); 


% Outside TRs
A=load(sprintf('mendelian_error_15x_%s_control_not_tr.csv',EVAL_THRESHOLD));
B=load(sprintf('mendelian_error_15x_%s_aou_not_tr.csv',EVAL_THRESHOLD));
C=load(sprintf('mendelian_error_30x_%s_aou_not_tr.csv',EVAL_THRESHOLD));
subplot(1,3,3); hold on;
[nrows,ncolumns]=size(A);
for i=[2:2:ncolumns]
	[nrows,ncolumns]=size(A); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=A(:,i)./(A(:,i)+A(:,i-1)); plot(X,Y,'.r');
	[nrows,ncolumns]=size(B); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=B(:,i)./(B(:,i)+B(:,i-1)); plot(X,Y,'.b');
	[nrows,ncolumns]=size(C); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=C(:,i)./(C(:,i)+C(:,i-1)); plot(X,Y,'ob');
endfor
ylabel('Mendelian error rate'); xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('Outside TRs'); grid on; axis([0,length(LABELS)+1,0,0.25]); axis square; set(gca,'fontsize',FONT_SIZE);
legend('15x control','15x AoU','30x AoU', 'location','southoutside'); 




% ------------------------------ De novo rate ----------------------------------
figure(3);

% All records
A=load(sprintf('denovo_15x_%s_control_all.csv',EVAL_THRESHOLD));
B=load(sprintf('denovo_15x_%s_aou_all.csv',EVAL_THRESHOLD));
C=load(sprintf('denovo_30x_%s_aou_all.csv',EVAL_THRESHOLD));
subplot(1,3,1); hold on;
[nrows,ncolumns]=size(A);
for i=[2:2:ncolumns]
	[nrows,ncolumns]=size(A); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=A(:,i-1)./A(:,i); plot(X,Y,'.r');
	[nrows,ncolumns]=size(B); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=B(:,i-1)./B(:,i); plot(X,Y,'.b');
	[nrows,ncolumns]=size(C); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=C(:,i-1)./C(:,i); plot(X,Y,'ob');
endfor
ylabel('De novo rate'); xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('All records'); grid on; axis([0,length(LABELS)+1,0,0.25]); axis square; set(gca,'fontsize',FONT_SIZE);
legend('15x control','15x AoU','30x AoU', 'location','southoutside'); 


% Inside TRs
A=load(sprintf('denovo_15x_%s_control_tr.csv',EVAL_THRESHOLD));
B=load(sprintf('denovo_15x_%s_aou_tr.csv',EVAL_THRESHOLD));
C=load(sprintf('denovo_30x_%s_aou_tr.csv',EVAL_THRESHOLD));
subplot(1,3,2); hold on;
[nrows,ncolumns]=size(A);
for i=[2:2:ncolumns]
	[nrows,ncolumns]=size(A); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=A(:,i-1)./A(:,i); plot(X,Y,'.r');
	[nrows,ncolumns]=size(B); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=B(:,i-1)./B(:,i); plot(X,Y,'.b');
	[nrows,ncolumns]=size(C); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=C(:,i-1)./C(:,i); plot(X,Y,'ob');
endfor
ylabel('De novo rate'); xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('Inside TRs'); grid on; axis([0,length(LABELS)+1,0,0.25]); axis square; set(gca,'fontsize',FONT_SIZE);
legend('15x control','15x AoU','30x AoU', 'location','southoutside'); 


% Outside TRs
A=load(sprintf('denovo_15x_%s_control_not_tr.csv',EVAL_THRESHOLD));
B=load(sprintf('denovo_15x_%s_aou_not_tr.csv',EVAL_THRESHOLD));
C=load(sprintf('denovo_30x_%s_aou_not_tr.csv',EVAL_THRESHOLD));
subplot(1,3,3); hold on;
[nrows,ncolumns]=size(A);
for i=[2:2:ncolumns]
	[nrows,ncolumns]=size(A); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=A(:,i-1)./A(:,i); plot(X,Y,'.r');
	[nrows,ncolumns]=size(B); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=B(:,i-1)./B(:,i); plot(X,Y,'.b');
	[nrows,ncolumns]=size(C); X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
	Y=C(:,i-1)./C(:,i); plot(X,Y,'ob');
endfor
ylabel('De novo rate'); xticks([1:length(MIN_N_SAMPLES)+1]); xticklabels(LABELS); xlabel('Min n. samples'); title('Outside TRs'); grid on; axis([0,length(LABELS)+1,0,0.25]); axis square; set(gca,'fontsize',FONT_SIZE);
legend('15x control','15x AoU','30x AoU', 'location','southoutside'); 
