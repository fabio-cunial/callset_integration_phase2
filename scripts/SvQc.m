
CALLER_ID=0;  % 0=pav, 1=pbsv, 2=sniffles

SV_TYPES_PRIMARY={'DEL','INS','DUP','INV'};
SV_TYPES_SECONDARY={'BND','SUB','UNK'};
MAX_COVERAGE=50;
MAX_NRECORDS=40000;
CALLER_NAMES={'PAV','PBSV','Sniffles'};


A=load('/Users/fcunial/Downloads/svqc/plot/counts.csv');
% Meaning of the columns of `counts.csv` (numbers ar offsets from the first 
% column, which contains the coverage):
%
% - >=20
%   - whole genome
%     - DEL
%       1 pav
%       2 pbsv
%       3 sniffles
%     - INS
%       4 pav
%       5 pbsv
%       6 sniffles
%     - DUP
%       7 pav
%       8 pbsv
%       9 sniffles
%     - INV
%       10 pav
%       11 pbsv
%       12 sniffles
%   - inside TR
%     - DEL
%       13 pav
%       14 pbsv
%       15 sniffles
%     - INS
%       16 pav
%       17 pbsv
%       18 sniffles
%     - DUP
%       19 pav
%       20 pbsv
%       21 sniffles
%     - INV
%       22 pav
%       23 pbsv
%       24 sniffles
%   - outside TR
%     - DEL
%       25 pav
%       26 pbsv
%       27 sniffles
%     - INS
%       28 pav
%       29 pbsv
%       30 sniffles
%     - DUP
%       31 pav
%       32 pbsv
%       33 sniffles
%     - INV
%       34 pav
%       35 pbsv
%       36 sniffles
%
% - >=50
%   - whole genome
%     - DEL
%       37 pav
%       38 pbsv
%       39 sniffles
%     - INS
%       40 pav
%       41 pbsv
%       42 sniffles
%     - DUP
%       43 pav
%       44 pbsv
%       45 sniffles
%     - INV
%       46 pav
%       47 pbsv
%       48 sniffles
%   - inside TR
%     - DEL
%       49 pav
%       50 pbsv
%       51 sniffles
%     - INS
%       52 pav
%       53 pbsv
%       54 sniffles
%     - DUP
%       55 pav
%       56 pbsv
%       57 sniffles
%     - INV
%       58 pav
%       59 pbsv
%       60 sniffles
%   - outside TR
%     - DEL
%       61 pav
%       62 pbsv
%       63 sniffles
%     - INS
%       64 pav
%       65 pbsv
%       66 sniffles
%     - DUP
%       67 pav
%       68 pbsv
%       69 sniffles
%     - INV
%       70 pav
%       71 pbsv
%       72 sniffles
%
% - Others
%   - whole genome
%    - BND
%       73 pav
%       74 pbsv
%       75 sniffles
%    - SUB
%       76 pav
%       77 pbsv
%       78 sniffles
%    - UNK
%       79 pav
%       80 pbsv
%       81 sniffles
%   - inside TR
%     - BND
%       82 pav
%       83 pbsv 
%       84 sniffles
%     - SUB
%       85 pav
%       86 pbsv
%       87 sniffles
%     - UNK
%       88 pav
%       89 pbsv
%       90 sniffles
%   - outside TR
%     - BND
%       91 pav
%       92 pbsv
%       93 sniffles
%     - SUB
%       94 pav  
%       95 pbsv
%       96 sniffles
%     - UNK
%       97 pav
%       98 pbsv
%       99 sniffles




figure(1);
x=A(:,1);


% >=20bp
subplot(3,3,1); hold on;
for i=[1,4,7,10]+CALLER_ID
    plot(x,A(:,1+i),'.');
end
xlabel('Coverage'); ylabel('Number of records'); title(sprintf('%s, >=20bp, whole genome', CALLER_NAMES{CALLER_ID+1})); axis([0 MAX_COVERAGE 0 MAX_NRECORDS]); axis square; grid on; set(gca,'fontsize',14); legend('DEL','INS','DUP','INV','location','eastoutside');
subplot(3,3,4); hold on;
for i=[13,16,19,22]+CALLER_ID
    plot(x,A(:,1+i),'.');
end
xlabel('Coverage'); ylabel('Number of records'); title(sprintf('%s, >=20bp, inside TR', CALLER_NAMES{CALLER_ID+1})); axis([0 MAX_COVERAGE 0 MAX_NRECORDS]); axis square; grid on; set(gca,'fontsize',14); legend('DEL','INS','DUP','INV','location','eastoutside');
subplot(3,3,7); hold on;
for i=[25,28,31,34]+CALLER_ID
    plot(x,A(:,1+i),'.');
end
xlabel('Coverage'); ylabel('Number of records'); title(sprintf('%s, >=20bp, outside TR', CALLER_NAMES{CALLER_ID+1})); axis([0 MAX_COVERAGE 0 MAX_NRECORDS]); axis square; grid on; set(gca,'fontsize',14); legend('DEL','INS','DUP','INV','location','eastoutside');


% >=50bp
subplot(3,3,2); hold on;
for i=[37,40,43,46]+CALLER_ID
    plot(x,A(:,1+i),'.');
end
xlabel('Coverage'); ylabel('Number of records'); title(sprintf('%s, >=50bp, whole genome', CALLER_NAMES{CALLER_ID+1})); axis([0 MAX_COVERAGE 0 MAX_NRECORDS]); axis square; grid on; set(gca,'fontsize',14); legend('DEL','INS','DUP','INV','location','eastoutside');
subplot(3,3,5); hold on;
for i=[49,52,55,58]+CALLER_ID
    plot(x,A(:,1+i),'.');
end
xlabel('Coverage'); ylabel('Number of records'); title(sprintf('%s, >=50bp, inside TR', CALLER_NAMES{CALLER_ID+1})); axis([0 MAX_COVERAGE 0 MAX_NRECORDS]); axis square; grid on; set(gca,'fontsize',14); legend('DEL','INS','DUP','INV','location','eastoutside');
subplot(3,3,8); hold on;
for i=[61,64,67,70]+CALLER_ID
    plot(x,A(:,1+i),'.');
end
xlabel('Coverage'); ylabel('Number of records'); title(sprintf('%s, >=50bp, outside TR', CALLER_NAMES{CALLER_ID+1})); axis([0 MAX_COVERAGE 0 MAX_NRECORDS]); axis square; grid on; set(gca,'fontsize',14); legend('DEL','INS','DUP','INV','location','eastoutside');


% Other types
subplot(3,3,3); hold on;
for i=[73,76,79]+CALLER_ID
    plot(x,A(:,1+i),'.');
end
xlabel('Coverage'); ylabel('Number of records'); title(sprintf('%s, Other types, whole genome', CALLER_NAMES{CALLER_ID+1})); axis([0 MAX_COVERAGE 0 MAX_NRECORDS]); axis square; grid on; set(gca,'fontsize',14); legend('BND','SUB','UNK','location','eastoutside');
subplot(3,3,6); hold on;
for i=[82,85,88]+CALLER_ID
    plot(x,A(:,1+i),'.');
end
xlabel('Coverage'); ylabel('Number of records'); title(sprintf('%s, Other types, inside TR', CALLER_NAMES{CALLER_ID+1})); axis([0 MAX_COVERAGE 0 MAX_NRECORDS]); axis square; grid on; set(gca,'fontsize',14); legend('BND','SUB','UNK','location','eastoutside');
subplot(3,3,9); hold on;
for i=[91,94,97]+CALLER_ID
    plot(x,A(:,1+i),'.');
end
xlabel('Coverage'); ylabel('Number of records'); title(sprintf('%s, Other types, outside TR', CALLER_NAMES{CALLER_ID+1})); axis([0 MAX_COVERAGE 0 MAX_NRECORDS]); axis square; grid on; set(gca,'fontsize',14); legend('BND','SUB','UNK','location','eastoutside');
