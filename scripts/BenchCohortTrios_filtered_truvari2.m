SOURCE_DIR='.';
FONT_SIZE=10;
DELTA=0.4;
LABELS={'V1','2','3','4','8','16','32','64','128','256','512','1024','2048'};
N_ARGS=length(LABELS);
N_TOTAL_ARGS=N_ARGS;
MAX_MENDELIAN_ERROR=0.2;
MAX_DE_NOVO=0.2;
PREFIX='filtered_truvari';


figure(1);




# --------------------------- Mendelian error ----------------------------------

# All, 07.
A=load(sprintf('%s/%s_all.csv',SOURCE_DIR,PREFIX));
[nrows,ncolumns]=size(A);

subplot(3,3,1); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
    plot(X,Y,'.');
endfor
title('All variants, 07.'); axis([0,N_ARGS+1,0,MAX_MENDELIAN_ERROR]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,4); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,N_TOTAL_ARGS*2+i);
    plot(X,Y,'.');
endfor
title('All variants, 07.'); axis([0,N_ARGS+1,0,MAX_DE_NOVO]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);


# TR, 07.
A=load(sprintf('%s/%s_tr.csv',SOURCE_DIR,PREFIX));
[nrows,ncolumns]=size(A);

subplot(3,3,2); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
    plot(X,Y,'.');
endfor
title('Inside TR, 07.'); axis([0,N_ARGS+1,0,MAX_MENDELIAN_ERROR]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,5); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,N_TOTAL_ARGS*2+i);
    plot(X,Y,'.');
endfor
title('Inside TR, 07.'); axis([0,N_ARGS+1,0,MAX_DE_NOVO]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);


# Not TR, 07.
A=load(sprintf('%s/%s_not_tr.csv',SOURCE_DIR,PREFIX));
[nrows,ncolumns]=size(A);

subplot(3,3,3); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
    plot(X,Y,'.');
endfor
title('Outside TR, 07.'); axis([0,N_ARGS+1,0,MAX_MENDELIAN_ERROR]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,6); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,N_TOTAL_ARGS*2+i);
    plot(X,Y,'.');
endfor
title('Outside TR, 07.'); axis([0,N_ARGS+1,0,MAX_DE_NOVO]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);





# -------------------------------- Counts --------------------------------------

A=load(sprintf('%s/%s_counts.csv',SOURCE_DIR,PREFIX));
[nrows,ncolumns]=size(A);

subplot(3,3,7); hold on;
plot([1:N_ARGS],A(:,2)./A(1,2),'.-');
title('All variants, 07.'); axis([0,N_ARGS+1,0,1]); grid on; ylabel('Fraction of calls in V1');
xticks([1:N_ARGS]); xticklabels(LABELS); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,8); hold on;
plot([1:N_ARGS],A(:,4)./A(1,4),'.-');
title('Inside TR, 07.'); axis([0,N_ARGS+1,0,1]); grid on; ylabel('Fraction of calls in V1');
xticks([1:N_ARGS]); xticklabels(LABELS); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,9); hold on;
plot([1:N_ARGS],A(:,6)./A(1,6),'.-');
title('Outside TR, 07.'); axis([0,N_ARGS+1,0,1]); grid on; ylabel('Fraction of calls in V1');
xticks([1:N_ARGS]); xticklabels(LABELS); set(gca,'fontsize',FONT_SIZE);
