SOURCE_DIR='.';
FONT_SIZE=14;
DELTA=0.4;
%LABELS={'V1','ND1','ND2', 'ND3','HS1','HS2', 'HS3','HS4','SZ1', 'SZ2','SS1','SS2', 'MQ1','MQ2','MQ3', 'MQ4'};
%LABELS={'V1','SZ1','SZ2','SS1','SS2','MQ1','MQ2','MQ3','MQ4'};
%LABELS={'V1','MAF_1','MAF_2','POA','WFA'};
%LABELS={'V1','MF_4','MF_6','PM_{20}','H_{20}','H','K_8'};
%LABELS={'V1','3-25-3','10-25-3'};
%LABELS={'V1','>=2','>=3','>=4'};
%LABELS={'TR_{V1}',    'aH','aTR','aTR_{<50}','aTR_{200}','aTR_{10k}','aTR_{>10k}','aTR_{>100}',    'SA','uSA','dSR_{50}','dSR_{150}','dSR_{>150}','hSR_{6}','hSR_{11}','hSR_{>11}','hSR_{>20}','iSR_{>10}','iSR_{>20}','qSR_{<50}','qSR_{<150}','qSR_{>150}','tSR_{<50}','tSR_{<150}','tSR_{>150}',    'SD','SD_{10k}',    'NU',    'G_{<15}','G_{20}','G_{25}','G_{30}','G_{55}','G_{60}','G_{65}','G_{70}','G_{75}','G_{80}','G_{85}','G_{>85}',    'VA'};
LABELS={'TR_{V1}',    'ML_1','ML_2','ML_4','ML_8','ML_{16}','ML_{32}','ML_{64}','ML_{128}',    'NR_1','NR_2','NR_4','NR_8','NR_{16}','NR_{32}','NR_{64}','NR_{128}',        'P_5','P_6','P_7','P_8','P_9','P_{95}','P_{96}','P_{97}','P_{98}','P_{99}','P_1',        'L_0','L_1','L_2','L_4','L_8','L_{16}','L_{32}',        'MF_1','MF_2','MF_3','MF_4','MF_5','MF_6','MF_7','MF_8','MF_9','MF_{95}','MF_{96}','MF_{97}','MF_{98}','MF_{99}','MF_{100}',        'SI_0','SI_1','SI_2','SI_4','SI_8','SI_{16}','SI_{32}',        'ST_0','ST_1','ST_2','ST_4','ST_8','ST_{16}','ST_{32}',        'IL_2','IL_4','IL_6','IL_8','IL_{10}','IL_{12}','IL_{14}','IL_{16}','IL_{18}','IL_{20}','IL_{40}','IL_{80}','IL_{160}','IL_{320}','IL_{640}','IL_{1280}','IL_{2560}' };
%DATASET='GIAB stratifications v3.6';
DATASET='Tandem Repeat Catalog v1.0';

N_ARGS=length(LABELS);
%N_TOTAL_ARGS=16;
N_TOTAL_ARGS=N_ARGS;
MAX_MENDELIAN_ERROR=0.35;
MAX_DE_NOVO=0.5;
PREFIX='tr_substrat';  %'filtered_truvari';  %'triohispanics';     %'triomerge';  %   'squish_trios';
TOTAL_RECORDS_07_50BP=1018505;  % All records, of any length: 2816177

% Columns:
% ${N_GOOD_ALT_i},${N_MERR_i},
% ...,
%
% ${DENOVO_i},
% ...


figure(1);

# TR, 07.
A=load(sprintf('%s/%s.csv',SOURCE_DIR,PREFIX));
[nrows,ncolumns]=size(A);

subplot(3,1,1); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,2*i)./(A(:,2*i)+A(:,2*i-1));
    plot(X,Y,'.');
endfor
title(sprintf('%s, 07.',DATASET)); axis([0,N_ARGS+1,0,MAX_MENDELIAN_ERROR]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('Mendelian error rate'); set(gca,'fontsize',FONT_SIZE);

subplot(3,1,2); hold on;
for i=[1:N_ARGS]
    X=ones(nrows,1).*i -DELTA/2 + rand(nrows,1).*DELTA;
    Y=A(:,N_TOTAL_ARGS*2+i);
    plot(X,Y,'.');
endfor
title(sprintf('%s, 07.',DATASET)); axis([0,N_ARGS+1,0,MAX_DE_NOVO]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
ylabel('De novo rate'); set(gca,'fontsize',FONT_SIZE);

subplot(3,1,3); hold on;
A=load(sprintf('%s_counts.txt',PREFIX));
bar(A./TOTAL_RECORDS_07_50BP);
title(sprintf('Fraction of all records >=50bp in 07 (total: %d)',TOTAL_RECORDS_07_50BP)); axis([0,N_ARGS+1,0,0.4]); grid on; 
xticks([1:N_ARGS]); xticklabels(LABELS);
set(gca,'fontsize',FONT_SIZE);
