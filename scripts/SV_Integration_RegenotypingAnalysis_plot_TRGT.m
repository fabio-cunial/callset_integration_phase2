FONT_SIZE=14;
LABELS={">=20bp",">=50bp"};
DELTA=0.4;





% ----------------------------- Precision/recall -------------------------------
figure(1);

% All records
subplot(1,3,1); hold on; 
A=load('./20bp/precision_recall_20bp_all.csv');
B=load('./50bp/precision_recall_50bp_all.csv');
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2);
    plot(X,P,'ob'); plot(X,R,'or');
    
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,3); R=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r');
    
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=B(i,1); R=B(i,2);
    plot(X,P,'ob'); plot(X,R,'or');
    
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=B(i,3); R=B(i,4);
    plot(X,P,'.b'); plot(X,R,'.r');
endfor
xticks([1:2]); xticklabels(LABELS); title('All records'); grid on; axis([0,length(LABELS)+1,0.8,1]); axis square; 
legend('15x precision, with TRGT','15x recall, with TRGT', '15x precision, no TRGT','15x recall, no TRGT', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);

% TR
subplot(1,3,2); hold on; 
A=load('./20bp/precision_recall_20bp_tr.csv');
B=load('./50bp/precision_recall_50bp_tr.csv');
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2);
    plot(X,P,'ob'); plot(X,R,'or');
    
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,3); R=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r');
    
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=B(i,1); R=B(i,2);
    plot(X,P,'ob'); plot(X,R,'or');
    
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=B(i,3); R=B(i,4);
    plot(X,P,'.b'); plot(X,R,'.r');
endfor
xticks([1:2]); xticklabels(LABELS); title('Inside TR'); grid on; axis([0,length(LABELS)+1,0.8,1]); axis square; 
legend('15x precision, with TRGT','15x recall, with TRGT', '15x precision, no TRGT','15x recall, no TRGT', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);

% Not TR
subplot(1,3,3); hold on; 
A=load('./20bp/precision_recall_20bp_not_tr.csv');
B=load('./50bp/precision_recall_50bp_not_tr.csv');
[nrows,ncolumns]=size(A);
for i=[1:nrows]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,1); R=A(i,2);
    plot(X,P,'ob'); plot(X,R,'or');
    
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=A(i,3); R=A(i,4);
    plot(X,P,'.b'); plot(X,R,'.r');
    
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=B(i,1); R=B(i,2);
    plot(X,P,'ob'); plot(X,R,'or');
    
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=B(i,3); R=B(i,4);
    plot(X,P,'.b'); plot(X,R,'.r');
endfor
xticks([1:2]); xticklabels(LABELS); title('Outside TR'); grid on; axis([0,length(LABELS)+1,0.8,1]); axis square; 
legend('15x precision, with TRGT','15x recall, with TRGT', '15x precision, no TRGT','15x recall, no TRGT', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);
