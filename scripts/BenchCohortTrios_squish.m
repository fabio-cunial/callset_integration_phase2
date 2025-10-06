SOURCE_DIR='.';
FONT_SIZE=16;
DELTA=0.4;
%LABELS={'V1','ND1','ND2', 'ND3','HS1','HS2', 'HS3','HS4','SZ1', 'SZ2','SS1','SS2', 'MQ1','MQ2','MQ3', 'MQ4'};
%LABELS={'V1','SZ1','SZ2','SS1','SS2','MQ1','MQ2','MQ3','MQ4'};
%LABELS={'V1','MAF_1','MAF_2','POA','WFA'};
%LABELS={'V1','MF_4','MF_6','PM_{20}','H_{20}','H','K_8'};
%LABELS={'V1','3-25-3','10-25-3'};
%LABELS={'V1','>=2','>=3','>=4'};
LABELS={'V1','H','TR','S','SD','SD_{10k}','NU','G_{<15}','G_{20}','G_{25}','G_{30}','G_{55}','G_{60}','G_{65}','G_{70}','G_{75}','G_{80}','G_{85}','G_{>85}'};
N_ARGS=length(LABELS);
%N_TOTAL_ARGS=16;
N_TOTAL_ARGS=19;
MAX_MENDELIAN_ERROR=0.2;
MAX_DE_NOVO=0.2;
PREFIX='tr_substrat';  %'filtered_truvari';  %'triohispanics';     %'triomerge';  %   'squish_trios';

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


# TR, 07.
A=load(sprintf('%s/%s_tr.csv',SOURCE_DIR,PREFIX));
[nrows,ncolumns]=size(A);

subplot(2,3,2); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
    plot(X,Y,'.');
endfor
title('Inside TR, 07.'); axis([0,N_ARGS+1,0,MAX_MENDELIAN_ERROR]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

subplot(2,3,5); hold on;
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

subplot(2,3,3); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
    plot(X,Y,'.');
endfor
title('Outside TR, 07.'); axis([0,N_ARGS+1,0,MAX_MENDELIAN_ERROR]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

subplot(2,3,6); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,N_TOTAL_ARGS*2+i);
    plot(X,Y,'.');
endfor
title('Outside TR, 07.'); axis([0,N_ARGS+1,0,MAX_DE_NOVO]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);
