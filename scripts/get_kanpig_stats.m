figure(1);

subplot(1,2,1); hold on;
A=load('workpackage1_table1.csv');
[yy,xx]=hist(A(:,3),50); 
bar(xx,yy,'facecolor','blue','edgecolor','blue');
A=load('workpackage9_table1.csv');
[yy,xx]=hist(A(:,3),50);
bar(xx,yy,'facecolor','red','edgecolor','red');
xlabel('% records marked as ALT'); ylabel('n. samples');
axis square; grid on; set(gca,'fontsize',14);
legend('Kanpig on intra-sample VCF','Kanpig on personalized cohort VCF');
title('Filtering induced by kanpig');


subplot(1,2,2); hold on;
A=load('workpackage1_table2.csv');
[yy,xx]=hist(A(:,3),50);
h=bar(xx,yy,'facecolor','blue','edgecolor','blue');
set(get(h, 'children'), 'facealpha', 0.5);
A=load('workpackage9_table2.csv');
[yy,xx]=hist(A(:,3),50);
h=bar(xx,yy,'facecolor','red','edgecolor','red');
set(get(h, 'children'), 'facealpha', 0.5);
xlabel('% records in autosomes marked as HET'); ylabel('n. samples'); title('Kanpig on intra-sample VCF');
axis square; grid on; set(gca,'fontsize',14);
%legend('Kanpig on intra-sample VCF','Kanpig on personalized cohort VCF');
title('HETs marked by kanpig');







% subplot(2,2,1);
% A=load('workpackage1_table1.csv');
% hist(A(:,3),50);
% xlabel('% records marked as ALT by kanpig'); ylabel('n. samples'); title('Kanpig on intra-sample VCF');
% axis square; grid on; set(gca,'fontsize',14);
%
% subplot(2,2,2);
% A=load('workpackage1_table2.csv');
% hist(A(:,3),50);
% xlabel('% records in autosomes marked as HET by kanpig'); ylabel('n. samples'); title('Kanpig on intra-sample VCF');
% axis square; grid on; set(gca,'fontsize',14);
%
% subplot(2,2,3);
% A=load('workpackage9_table1.csv');
% hist(A(:,3),50);
% xlabel('% records marked as ALT by kanpig'); ylabel('n. samples'); title('Kanpig on personalized cohort VCF');
% axis square; grid on; set(gca,'fontsize',14);
%
% subplot(2,2,4);
% A=load('workpackage9_table2.csv');
% hist(A(:,3),50);
% xlabel('% records in autosomes marked as HET by kanpig'); ylabel('n. samples'); title('Kanpig on personalized cohort VCF');
% axis square; grid on; set(gca,'fontsize',14);
