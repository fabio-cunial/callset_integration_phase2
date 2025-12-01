FONT_SIZE=14;
LABELS={">=20bp",">=50bp"};
DELTA=0.4;





% ----------------------------- Precision/recall -------------------------------
figure(1);

% All records
subplot(1,3,1); hold on; 
AH=load('./regenotype_20/precision_recall_15x_20bp_haps_all.csv');
AH_90=load('./regenotype_20_collapsed_90/precision_recall_15x_20bp_haps_all.csv');
AR=load('./regenotype_20/precision_recall_15x_20bp_records_all.csv');
BH=load('./regenotype_50/precision_recall_15x_50bp_haps_all.csv');
BH_90=load('./regenotype_50_collapsed_90/precision_recall_15x_50bp_haps_all.csv');
BR=load('./regenotype_50/precision_recall_15x_50bp_records_all.csv');
[nrows,ncolumns]=size(AH);
for i=[1:nrows]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=AH(i,1); R=AH(i,2);
    plot(X,P,'ob'); plot(X,R,'or');
    
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=AH_90(i,1); R=AH_90(i,2);
    plot(X,P,'sb'); plot(X,R,'sr');
    
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=AR(i,1); R=AR(i,2);
    plot(X,P,'.b'); plot(X,R,'.r');
    
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=BH(i,1); R=BH(i,2);
    plot(X,P,'ob'); plot(X,R,'or');
    
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=BH_90(i,1); R=BH_90(i,2);
    plot(X,P,'sb'); plot(X,R,'sr');
    
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=BR(i,1); R=BR(i,2);
    plot(X,P,'.b'); plot(X,R,'.r');
endfor
xticks([1:2]); xticklabels(LABELS); title('All records'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; 
legend('15x precision, haps','15x recall, haps', '15x precision, haps collapsed 90%','15x recall, haps collapsed 90%', '15x precision, records','15x recall, records', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);


% TR
subplot(1,3,2); hold on; 
AR=load('./regenotype_20/precision_recall_15x_20bp_records_tr.csv');
BR=load('./regenotype_50/precision_recall_15x_50bp_records_tr.csv');
[nrows,ncolumns]=size(AR);
for i=[1:nrows]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=AR(i,1); R=AR(i,2);
    plot(X,P,'.b'); plot(X,R,'.r');
    
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=BR(i,1); R=BR(i,2);
    plot(X,P,'.b'); plot(X,R,'.r');
endfor
xticks([1:2]); xticklabels(LABELS); title('Inside TR'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; 
legend('15x precision, records','15x recall, records', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);


% NOT TR
subplot(1,3,3); hold on; 
AR=load('./regenotype_20/precision_recall_15x_20bp_records_not_tr.csv');
BR=load('./regenotype_50/precision_recall_15x_50bp_records_not_tr.csv');
[nrows,ncolumns]=size(AR);
for i=[1:nrows]
    X=1 -DELTA/2 + rand(1,1).*DELTA;
    P=AR(i,1); R=AR(i,2);
    plot(X,P,'.b'); plot(X,R,'.r');
    
    X=2 -DELTA/2 + rand(1,1).*DELTA;
    P=BR(i,1); R=BR(i,2);
    plot(X,P,'.b'); plot(X,R,'.r');
endfor
xticks([1:2]); xticklabels(LABELS); title('Outside TR'); grid on; axis([0,length(LABELS)+1,0,1]); axis square; 
legend('15x precision, records','15x recall, records', 'location','southoutside'); set(gca,'fontsize',FONT_SIZE);



% ----------------------------- Mendelian error --------------------------------
figure(2);


% All records
subplot(1,3,1); hold on;
AH=load('./regenotype_20/mendelian_error_15x_20bp_aou_haps_all.csv');
AH_90=load('./regenotype_20_collapsed_90/mendelian_error_15x_20bp_aou_haps_all.csv');
AR=load('./regenotype_20/mendelian_error_15x_20bp_aou_records_all.csv');
BH=load('./regenotype_50/mendelian_error_15x_50bp_aou_haps_all.csv');
BH_90=load('./regenotype_50_collapsed_90/mendelian_error_15x_50bp_aou_haps_all.csv');
BR=load('./regenotype_50/mendelian_error_15x_50bp_aou_records_all.csv');
[nrows,ncolumns]=size(AH);
for i=[2:2:ncolumns]
    [nrows,ncolumns]=size(AH);
    X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=AH(:,i)./(AH(:,i)+AH(:,i-1)); plot(X,Y,'ob');
    
    X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=AH_90(:,i)./(AH_90(:,i)+AH_90(:,i-1)); plot(X,Y,'sb');
    
    X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=AR(:,i)./(AR(:,i)+AR(:,i-1)); plot(X,Y,'.b');

    X=2*ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=BH(:,i)./(BH(:,i)+BH(:,i-1)); plot(X,Y,'ob');
    
    X=2*ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=BH_90(:,i)./(BH_90(:,i)+BH_90(:,i-1)); plot(X,Y,'sb');
    
    X=2*ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=BR(:,i)./(BR(:,i)+BR(:,i-1)); plot(X,Y,'.b');
endfor
ylabel('Mendelian error rate'); xticks([1:2]); xticklabels(LABELS); title('All records'); grid on; axis([0,length(LABELS)+1,0,0.25]); axis square; set(gca,'fontsize',FONT_SIZE);
legend('15x AoU, haps', '15x AoU, haps collapsed 90%', '15x AoU, records', 'location','southoutside');


% TR
subplot(1,3,2); hold on;
AR=load('./regenotype_20/mendelian_error_15x_20bp_aou_records_tr.csv');
BR=load('./regenotype_50/mendelian_error_15x_50bp_aou_records_tr.csv');
[nrows,ncolumns]=size(AR);
for i=[2:2:ncolumns]
    [nrows,ncolumns]=size(AR);
    X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=AR(:,i)./(AR(:,i)+AR(:,i-1)); plot(X,Y,'.b');

    X=2*ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=BR(:,i)./(BR(:,i)+BR(:,i-1)); plot(X,Y,'.b');
endfor
ylabel('Mendelian error rate'); xticks([1:2]); xticklabels(LABELS); title('Inside TR'); grid on; axis([0,length(LABELS)+1,0,0.25]); axis square; set(gca,'fontsize',FONT_SIZE);
legend('15x AoU, records', 'location','southoutside');


% NOT TR
subplot(1,3,3); hold on;
AR=load('./regenotype_20/mendelian_error_15x_20bp_aou_records_not_tr.csv');
BR=load('./regenotype_50/mendelian_error_15x_50bp_aou_records_not_tr.csv');
[nrows,ncolumns]=size(AR);
for i=[2:2:ncolumns]
    [nrows,ncolumns]=size(AR);
    X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=AR(:,i)./(AR(:,i)+AR(:,i-1)); plot(X,Y,'.b');

    X=2*ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=BR(:,i)./(BR(:,i)+BR(:,i-1)); plot(X,Y,'.b');
endfor
ylabel('Mendelian error rate'); xticks([1:2]); xticklabels(LABELS); title('Outside TR'); grid on; axis([0,length(LABELS)+1,0,0.25]); axis square; set(gca,'fontsize',FONT_SIZE);
legend('15x AoU, records', 'location','southoutside');





% ------------------------------ De novo rate ----------------------------------
figure(3);

% All records
subplot(1,3,1); hold on;
AH=load('./regenotype_20/denovo_15x_20bp_aou_haps_all.csv');
AH_90=load('./regenotype_20_collapsed_90/denovo_15x_20bp_aou_haps_all.csv');
AR=load('./regenotype_20/denovo_15x_20bp_aou_records_all.csv');
BH=load('./regenotype_50/denovo_15x_50bp_aou_haps_all.csv');
BH_90=load('./regenotype_50_collapsed_90/denovo_15x_50bp_aou_haps_all.csv');
BR=load('./regenotype_50/denovo_15x_50bp_aou_records_all.csv');
[nrows,ncolumns]=size(AH);
for i=[2:2:ncolumns]
    [nrows,ncolumns]=size(AH);
    X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=AH(:,i-1)./AH(:,i); plot(X,Y,'ob');
    
    X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=AH_90(:,i-1)./AH_90(:,i); plot(X,Y,'sb');
    
    X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=AR(:,i-1)./AR(:,i); plot(X,Y,'.b');

    X=2*ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=BH(:,i-1)./BH(:,i); plot(X,Y,'ob');
    
    X=2*ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=BH_90(:,i-1)./BH_90(:,i); plot(X,Y,'sb');
    
    X=2*ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=BR(:,i-1)./BR(:,i); plot(X,Y,'.b');
endfor
ylabel('De novo rate'); xticks([1:2]); xticklabels(LABELS); title('All records'); grid on; axis([0,length(LABELS)+1,0,0.25]); axis square; set(gca,'fontsize',FONT_SIZE);
legend('15x AoU, haps', '15x AoU, haps collapsed 90%', '15x AoU, records', 'location','southoutside');


% TR
subplot(1,3,2); hold on;
AR=load('./regenotype_20/denovo_15x_20bp_aou_records_tr.csv');
BR=load('./regenotype_50/denovo_15x_50bp_aou_records_tr.csv');
[nrows,ncolumns]=size(AH);
for i=[2:2:ncolumns]
    [nrows,ncolumns]=size(AR);
    
    X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=AR(:,i-1)./AR(:,i); plot(X,Y,'.b');

    X=2*ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=BR(:,i-1)./BR(:,i); plot(X,Y,'.b');
endfor
ylabel('De novo rate'); xticks([1:2]); xticklabels(LABELS); title('Inside TR'); grid on; axis([0,length(LABELS)+1,0,0.25]); axis square; set(gca,'fontsize',FONT_SIZE);
legend('15x AoU, records', 'location','southoutside');


% NOT TR
subplot(1,3,3); hold on;
AR=load('./regenotype_20/denovo_15x_20bp_aou_records_not_tr.csv');
BR=load('./regenotype_50/denovo_15x_50bp_aou_records_not_tr.csv');
[nrows,ncolumns]=size(AH);
for i=[2:2:ncolumns]
    [nrows,ncolumns]=size(AR);
    
    X=ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=AR(:,i-1)./AR(:,i); plot(X,Y,'.b');

    X=2*ones(nrows,1).*i/2 -DELTA/2 + rand(nrows,1).*DELTA;
    Y=BR(:,i-1)./BR(:,i); plot(X,Y,'.b');
endfor
ylabel('De novo rate'); xticks([1:2]); xticklabels(LABELS); title('Outside TR'); grid on; axis([0,length(LABELS)+1,0,0.25]); axis square; set(gca,'fontsize',FONT_SIZE);
legend('15x AoU, records', 'location','southoutside');
