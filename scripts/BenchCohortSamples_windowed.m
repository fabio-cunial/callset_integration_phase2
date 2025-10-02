SOURCE_DIR='.';
FONT_SIZE=18;
MAX_WINDOW=98;
WINDOW_SIZE=2500000;
CHROMOSOME_LENGTH=250000000;
X=[1:MAX_WINDOW+1]*WINDOW_SIZE;


kanpig=load('kanpig.csv'); 
f07=load('07.csv');
f09=load('09.csv');
cohort_merged_07=load('cohort_merged_07.csv');
cohort_merged_09=load('cohort_merged_09.csv');
cohort_regenotyped_07=load('cohort_regenotyped_07.csv');
cohort_regenotyped_09=load('cohort_regenotyped_09.csv');


figure(1);


subplot(4,2,1); hold on;
bar(X,kanpig(:,1));
title('Precision of single-sample kanpig'); axis([0,CHROMOSOME_LENGTH,0,1]); grid on; xlabel('chr1 POS'); ylabel('P_{kanpig}'); set(gca,'fontsize',FONT_SIZE);
subplot(4,2,2); hold on;
bar(X,kanpig(:,2));
title('Recall of single-sample kanpig'); axis([0,CHROMOSOME_LENGTH,0,1]); grid on; xlabel('chr1 POS'); ylabel('R_{kanpig}'); set(gca,'fontsize',FONT_SIZE);



subplot(4,2,3); hold on;
plot(X,f07(:,1)-kanpig(:,1),'.');
plot(X,f09(:,1)-kanpig(:,1),'o');
title('Precision \Delta of XGBoost'); axis([0,CHROMOSOME_LENGTH,-0.5,1]); grid on; xlabel('chr1 POS'); ylabel('P-P_{kanpig}'); set(gca,'fontsize',FONT_SIZE);
legend('07','09','location','northeastoutside');
subplot(4,2,4); hold on;
plot(X,f07(:,2)-kanpig(:,2),'.');
plot(X,f09(:,2)-kanpig(:,2),'o');
title('Recall \Delta of XGBoost'); axis([0,CHROMOSOME_LENGTH,-0.5,1]); grid on; xlabel('chr1 POS'); ylabel('R-R_{kanpig}'); set(gca,'fontsize',FONT_SIZE);
legend('07','09','location','northeastoutside');



subplot(4,2,5); hold on;
plot(X,cohort_merged_07(:,1)-f07(:,1),'.');
plot(X,cohort_merged_09(:,1)-f09(:,1),'o');
title('Precision \Delta of truvari collapse'); axis([0,CHROMOSOME_LENGTH,-0.5,0.5]); grid on; xlabel('chr1 POS'); ylabel('P-P_{xgboost}'); set(gca,'fontsize',FONT_SIZE);
legend('07','09','location','northeastoutside');
subplot(4,2,6); hold on;
plot(X,cohort_merged_07(:,2)-f07(:,2),'.');
plot(X,cohort_merged_09(:,2)-f09(:,2),'o');
title('Recall \Delta of truvari collapse'); axis([0,CHROMOSOME_LENGTH,-0.5,0.5]); grid on; xlabel('chr1 POS'); ylabel('R-R_{xgboost}'); set(gca,'fontsize',FONT_SIZE);
legend('07','09','location','northeastoutside');




subplot(4,2,7); hold on;
plot(X,cohort_regenotyped_07(:,1)-cohort_merged_07(:,1),'.');
plot(X,cohort_regenotyped_09(:,1)-cohort_merged_09(:,1),'o');
title('Precision \Delta of re-genotyping'); axis([0,CHROMOSOME_LENGTH,-0.8,0.6]); grid on; xlabel('chr1 POS'); ylabel('P-P_{truvari}'); set(gca,'fontsize',FONT_SIZE);
legend('07','09','location','northeastoutside');
subplot(4,2,8); hold on;
plot(X,cohort_regenotyped_07(:,2)-cohort_merged_07(:,2),'.');
plot(X,cohort_regenotyped_09(:,2)-cohort_merged_09(:,2),'o');
title('Recall \Delta of re-genotyping'); axis([0,CHROMOSOME_LENGTH,-0.8,0.6]); grid on; xlabel('chr1 POS'); ylabel('R-R_{truvari}'); set(gca,'fontsize',FONT_SIZE);
legend('07','09','location','northeastoutside');





%
% subplot(1,2,1); hold on;
% plot(kanpig(:,1),'.');
% plot(f07(:,1)),'.';
% plot(f09(:,1),'.');
% plot(cohort_merged_07(:,1),'.');
% plot(cohort_merged_09(:,1),'.');
% plot(cohort_regenotyped_07(:,1),'.');
% plot(cohort_regenotyped_09(:,1),'.');
% title('Precision - chr1'); axis([0,MAX_WINDOW,0,1]); grid on; xlabel('chr1 POS'); ylabel('Precision'); set(gca,'fontsize',FONT_SIZE);
%
% subplot(1,2,2); hold on;
% plot(kanpig(:,2),'.');
% plot(f07(:,2),'.');
% plot(f09(:,2),'.');
% plot(cohort_merged_07(:,2),'.');
% plot(cohort_merged_09(:,2),'.');
% plot(cohort_regenotyped_07(:,2),'.');
% plot(cohort_regenotyped_09(:,2),'.');
% title('Recall - chr1'); axis([0,MAX_WINDOW,0,1]); grid on; xlabel('chr1 POS'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);
