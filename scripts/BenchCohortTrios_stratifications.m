SOURCE_DIR='.';
FONT_SIZE=16;
DELTA=0.4;
LABELS={'V1','D','N'};
N_ARGS=length(LABELS);
N_TOTAL_ARGS=3;
MAX_MENDELIAN_ERROR=0.12;
MAX_DE_NOVO=0.12;
PREFIX='stratifications';

% Columns:
% ${N_GOOD_ALT_i},${N_MERR_i},
% ...,
%
% ${DENOVO_i},
% ...


figure(1);

# All, 07.
A=load(sprintf('%s/%s_all.csv',SOURCE_DIR,PREFIX));

[nrows,ncolumns]=size(A);

subplot(2,3,1); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
    plot(X,Y,'.');
endfor
title('All variants, 07.'); axis([0,N_ARGS+1,0,MAX_MENDELIAN_ERROR]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

subplot(2,3,4); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,N_TOTAL_ARGS*2+i);
    plot(X,Y,'.');
endfor
title('All variants, 07.'); axis([0,N_ARGS+1,0,MAX_DE_NOVO]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);
