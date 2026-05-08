LABELS={'N','T','TI','TC','I','P'};




figure(1);

subplot(2,2,1); hold on;
A=load('remap_classification_0_10000000.csv');
B=sum(A);
bar(100*B(2:end)./B(1,1)); 
xticks(1:6); xticklabels(LABELS); ylabel('% records'); xlabel('Class'); title('INS / Remap classification / All'); grid on; set(gca,'fontsize',14);
ylim([0 100]);

subplot(2,2,2); hold on;
A=load('remap_classification_10000_20000.csv');
B=sum(A);
bar(100*B(2:end)./B(1,1)); 
xticks(1:6); xticklabels(LABELS); ylabel('% records'); xlabel('Class'); title('INS / Remap classification / 10kb..20kb'); grid on; set(gca,'fontsize',14);
ylim([0 100]);

subplot(2,2,3); hold on;
A=load('remap_classification_20000_50000.csv');
B=sum(A);
bar(100*B(2:end)./B(1,1)); 
xticks(1:6); xticklabels(LABELS); ylabel('% records'); xlabel('Class'); title('INS / Remap classification / 20kb..50kb'); grid on; set(gca,'fontsize',14);
ylim([0 100]);

subplot(2,2,4); hold on;
A=load('remap_classification_50000_100000.csv');
B=sum(A);
bar(100*B(2:end)./B(1,1)); 
xticks(1:6); xticklabels(LABELS); ylabel('% records'); xlabel('Class'); title('INS / Remap classification / 50kb..100kb'); grid on; set(gca,'fontsize',14);
ylim([0 100]);




figure(2);

subplot(2,2,1); hold on;
A=load('remap_sniffles.csv');
B=sum(A);
bar(100*B(2:end)./B(1,1)); 
xticks(1:6); xticklabels(LABELS); ylabel('% of all sniffles records'); xlabel('Class'); title('Sniffles'); grid on; set(gca,'fontsize',14);
ylim([0 100]);

subplot(2,2,2); hold on;
A=load('remap_pbsv.csv');
B=sum(A);
bar(100*B(2:end)./B(1,1)); 
xticks(1:6); xticklabels(LABELS); ylabel('% of all pbsv records'); xlabel('Class'); title('PBSV'); grid on; set(gca,'fontsize',14);
ylim([0 100]);

subplot(2,2,3); hold on;
A=load('remap_pav.csv');
B=sum(A);
bar(100*B(2:end)./B(1,1)); 
xticks(1:6); xticklabels(LABELS); ylabel('% of all pav records'); xlabel('Class'); title('PAV'); grid on; set(gca,'fontsize',14);
ylim([0 100]);




figure(3);

subplot(1,3,1); hold on;
A=load('remap_tandems.csv');
B=sum(A);
bar(100*B(2:end)./B(1,1)); 
xticks(1:3); xticklabels({'Sniffles','PBSV','PAV'}); ylabel('% of all tandem records'); xlabel('Caller'); title('Tandem records by caller'); grid on; set(gca,'fontsize',14);
ylim([0 100]);

subplot(1,3,2); hold on;
A=load('remap_tandem_lengths.csv');
B=sum(A);
bar(100*B(2:end)./B(1,1)); 
xticks(1:4); xticklabels({'10-20kb','20-50kb','50-100kb','>100kb'}); ylabel('% of all tandem records'); xlabel('Length'); title('Length distribution of tandem records'); grid on; set(gca,'fontsize',14);
ylim([0 100]);

subplot(1,3,3); hold on;
A=load('remap_tandem_by_length.csv');
B=sum(A);
bar([100*B(2)./B(1), 100*B(4)./B(3), 100*B(6)./B(5), 100*B(8)./B(7)]); 
xticks(1:4); xticklabels({'10-20kb','20-50kb','50-100kb','>100kb'}); ylabel('% of records marked as tandem'); xlabel('Length'); title('Fraction of tandem records at each length'); grid on; set(gca,'fontsize',14);
ylim([0 100]);