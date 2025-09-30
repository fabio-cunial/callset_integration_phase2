INPUT_DIR="/Users/fcunial/Downloads/denovo_by_region";
SAMPLES={'1004931','1402381','1665275','1806012','3498199'};
DELTA=0.4;
FONT_SIZE=14;
MAX_MENDELIAN_ERROR=0.15;


subplot(1,3,1); hold on;
for i=[1:length(SAMPLES)]
    X=[1:23]-DELTA/2+rand(1,23).*DELTA;
    Y=zeros(1,23);
    for j=[1:22]
        A=load(sprintf('%s/%s_all_chr%d.tsv',INPUT_DIR,SAMPLES{i},j));
        Y(j)=sum(A(:,4))/sum(A(:,3));
    endfor
    A=load(sprintf('%s/%s_all_chrX.tsv',INPUT_DIR,SAMPLES{i}));
    Y(23)=sum(A(:,4))/sum(A(:,3));
    plot(X,Y,'.');
endfor
title('All calls'); axis([0,23,0,MAX_MENDELIAN_ERROR]); grid on; set(gca,'fontsize',FONT_SIZE);
xticks([1:24]); xlabel('Chromosome ID'); ylabel('Mendelian error rate');


subplot(1,3,2); hold on;
for i=[1:length(SAMPLES)]
    X=[1:23]-DELTA/2+rand(1,23).*DELTA;
    Y=zeros(1,23);
    for j=[1:22]
        A=load(sprintf('%s/%s_tr_chr%d.tsv',INPUT_DIR,SAMPLES{i},j));
        Y(j)=sum(A(:,4))/sum(A(:,3));
    endfor
    A=load(sprintf('%s/%s_tr_chrX.tsv',INPUT_DIR,SAMPLES{i}));
    Y(23)=sum(A(:,4))/sum(A(:,3));
    plot(X,Y,'.');
endfor
title('Inside TR'); axis([0,23,0,MAX_MENDELIAN_ERROR]); grid on; set(gca,'fontsize',FONT_SIZE);
xticks([1:24]); xlabel('Chromosome ID'); ylabel('Mendelian error rate');


subplot(1,3,3); hold on;
for i=[1:length(SAMPLES)]
    X=[1:23]-DELTA/2+rand(1,23).*DELTA;
    Y=zeros(1,23);
    for j=[1:22]
        A=load(sprintf('%s/%s_not_tr_chr%d.tsv',INPUT_DIR,SAMPLES{i},j));
        Y(j)=sum(A(:,4))/sum(A(:,3));
    endfor
    A=load(sprintf('%s/%s_not_tr_chrX.tsv',INPUT_DIR,SAMPLES{i}));
    Y(23)=sum(A(:,4))/sum(A(:,3));
    plot(X,Y,'.');
endfor
title('Outside TR'); axis([0,23,0,MAX_MENDELIAN_ERROR]); grid on; set(gca,'fontsize',FONT_SIZE);
xticks([1:24]); xlabel('Chromosome ID'); ylabel('Mendelian error rate');


% -----------------------

figure(2); hold on;
for i=[1:length(SAMPLES)]
    A=load(sprintf('%s/%s_all_chr1.tsv',INPUT_DIR,SAMPLES{i}));
    plot((A(:,2)+A(:,1))/2,A(:,4)./A(:,3),'.');
endfor
ylabel('De novo rate per kanpig window'); xlabel('Center of a kanpig window (chr1)'); grid on; set(gca,'fontsize',FONT_SIZE);




figure(3);
QUANTUM=2500000;  % Arbitrary
for i=[1:3]
    subplot(5,1,i); hold on;
    B=zeros(100,3);
    A=load(sprintf('%s/%s_all_chr1.tsv',INPUT_DIR,SAMPLES{i}));
    for k=[1:length(A)]
        j=floor( (A(k,2)+A(k,1))/(2*QUANTUM) )+1;
        if (j<=100) 
            B(j,1)=B(j,1)+A(k,3);
            B(j,2)=B(j,2)+A(k,4);
            B(j,3)=B(j,3)+A(k,5);
        endif
    endfor
    bar([1:100].*QUANTUM,B(:,2)./B(:,1),'stacked');
    xlabel('chr1 POS'); title(sprintf('Mendelian error rate, child %s',SAMPLES{i})); grid on; axis([0,250000000,0,0.4]); set(gca,'fontsize',FONT_SIZE);
endfor


# Computed with bcftools view on every chunk of length 25000 of
# `truvari_collapsed_for_kanpig.vcf.gz` (stringent).
ALL_CALLS_07_CHR1=load(sprintf('%s/all_calls_07_chr1.txt',INPUT_DIR));
subplot(5,1,4);
bar([1:100].*QUANTUM,ALL_CALLS_07_CHR1);
xlabel('chr1 POS'); title('Total calls in the stringent inter-sample VCF'); grid on; axis([0,250000000,0,12000]); set(gca,'fontsize',FONT_SIZE);

ALL_CALLS_09_CHR1=load(sprintf('%s/all_calls_09_chr1.txt',INPUT_DIR));
subplot(5,1,5);
bar([1:100].*QUANTUM,ALL_CALLS_09_CHR1);
xlabel('chr1 POS'); title('Total calls in the lenient inter-sample VCF'); grid on; axis([0,250000000,0,25000]); set(gca,'fontsize',FONT_SIZE);




%
% centers=zeros(23,2);
% radii=zeros(23,1);
% for i=[1:22]
%     A=load(sprintf('%s/%s_tr_chr%d.tsv',INPUT_DIR,PREFIX,i));
%     centers(i,1)=sum(A(:,4))/sum(A(:,3));
%     radii(i)=sum(A(:,3));
%     A=load(sprintf('%s/%s_not_tr_chr%d.tsv',INPUT_DIR,PREFIX,i));
%     centers(i,2)=sum(A(:,4))/sum(A(:,3));
%     radii(i)=radii(i)+sum(A(:,3));
% endfor
% A=load(sprintf('%s/%s_tr_chrX.tsv',INPUT_DIR,PREFIX));
% centers(23,1)=sum(A(:,4))/sum(A(:,3));
% radii(23)=sum(A(:,3));
% A=load(sprintf('%s/%s_not_tr_chrX.tsv',INPUT_DIR,PREFIX));
% centers(23,2)=sum(A(:,4))/sum(A(:,3));
% radii(23)=radii(i)+sum(A(:,3));
% plot(centers(:,1),centers(:,2),'o');
