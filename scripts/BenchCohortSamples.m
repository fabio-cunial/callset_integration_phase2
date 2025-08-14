SOURCE_DIR='.';
FONT_SIZE=18;


# All calls
kanpig=load('kanpig_all.csv'); 
f07=load('07_all.csv');
f09=load('09_all.csv');
cohort_merged_07=load('cohort_merged_07_all.csv');
cohort_merged_09=load('cohort_merged_09_all.csv');
cohort_regenotyped_07=load('cohort_regenotyped_07_all.csv');
cohort_regenotyped_09=load('cohort_regenotyped_09_all.csv');
subplot(3,2,1); hold on;
P07=[kanpig(:,1),f07(:,1),cohort_merged_07(:,1),cohort_regenotyped_07(:,1)];
R07=[kanpig(:,2),f07(:,2),cohort_merged_07(:,2),cohort_regenotyped_07(:,2)];
plot(P07(:,1),R07(:,1),'.r');
plot(P07(:,2),R07(:,2),'.g');
plot(P07(:,3),R07(:,3),'.b');
plot(P07(:,4),R07(:,4),'.y');
title('All calls, 07'); axis([0.8,1,0.4,0.8]); axis square; grid on; xlabel('Precision'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);
subplot(3,2,2); hold on;
P09=[kanpig(:,1),f09(:,1),cohort_merged_09(:,1),cohort_regenotyped_09(:,1)];
R09=[kanpig(:,2),f09(:,2),cohort_merged_09(:,2),cohort_regenotyped_09(:,2)];
plot(P09(:,1),R09(:,1),'.r');
plot(P09(:,2),R09(:,2),'.g');
plot(P09(:,3),R09(:,3),'.b');
plot(P09(:,4),R09(:,4),'.y');
title('All calls, 09'); axis([0.8,1,0.4,0.8]); axis square; grid on; xlabel('Precision'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);


# Inside TRs
kanpig=load('kanpig_tr.csv'); 
f07=load('07_tr.csv');
f09=load('09_tr.csv');
cohort_merged_07=load('cohort_merged_07_tr.csv');
cohort_merged_09=load('cohort_merged_09_tr.csv');
cohort_regenotyped_07=load('cohort_regenotyped_07_tr.csv');
cohort_regenotyped_09=load('cohort_regenotyped_09_tr.csv');
subplot(3,2,3); hold on;
P07=[kanpig(:,1),f07(:,1),cohort_merged_07(:,1),cohort_regenotyped_07(:,1)];
R07=[kanpig(:,2),f07(:,2),cohort_merged_07(:,2),cohort_regenotyped_07(:,2)];
plot(P07(:,1),R07(:,1),'.r');
plot(P07(:,2),R07(:,2),'.g');
plot(P07(:,3),R07(:,3),'.b');
plot(P07(:,4),R07(:,4),'.y');
title('Inside TR, 07'); axis([0.8,1,0.4,0.8]); axis square; grid on; xlabel('Precision'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);
subplot(3,2,4); hold on;
P09=[kanpig(:,1),f09(:,1),cohort_merged_09(:,1),cohort_regenotyped_09(:,1)];
R09=[kanpig(:,2),f09(:,2),cohort_merged_09(:,2),cohort_regenotyped_09(:,2)];
plot(P09(:,1),R09(:,1),'.r');
plot(P09(:,2),R09(:,2),'.g');
plot(P09(:,3),R09(:,3),'.b');
plot(P09(:,4),R09(:,4),'.y');
title('Inside TR, 09'); axis([0.8,1,0.4,0.8]); axis square; grid on; xlabel('Precision'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);


# Outside TRs
kanpig=load('kanpig_not_tr.csv'); 
f07=load('07_not_tr.csv');
f09=load('09_not_tr.csv');
cohort_merged_07=load('cohort_merged_07_not_tr.csv');
cohort_merged_09=load('cohort_merged_09_not_tr.csv');
cohort_regenotyped_07=load('cohort_regenotyped_07_not_tr.csv');
cohort_regenotyped_09=load('cohort_regenotyped_09_not_tr.csv');
subplot(3,2,5); hold on;
P07=[kanpig(:,1),f07(:,1),cohort_merged_07(:,1),cohort_regenotyped_07(:,1)];
R07=[kanpig(:,2),f07(:,2),cohort_merged_07(:,2),cohort_regenotyped_07(:,2)];
plot(P07(:,1),R07(:,1),'.r');
plot(P07(:,2),R07(:,2),'.g');
plot(P07(:,3),R07(:,3),'.b');
plot(P07(:,4),R07(:,4),'.y');
title('Outside TR, 07'); axis([0.8,1,0.6,1]); axis square; grid on; xlabel('Precision'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);
subplot(3,2,6); hold on;
P09=[kanpig(:,1),f09(:,1),cohort_merged_09(:,1),cohort_regenotyped_09(:,1)];
R09=[kanpig(:,2),f09(:,2),cohort_merged_09(:,2),cohort_regenotyped_09(:,2)];
plot(P09(:,1),R09(:,1),'.r');
plot(P09(:,2),R09(:,2),'.g');
plot(P09(:,3),R09(:,3),'.b');
plot(P09(:,4),R09(:,4),'.y');
title('Outside TR, 09'); axis([0.8,1,0.6,1]); axis square; grid on; xlabel('Precision'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);
