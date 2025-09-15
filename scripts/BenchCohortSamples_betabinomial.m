SOURCE_DIR='.';
FONT_SIZE=18;


# All calls
kanpig=load('kanpig_all.csv'); 
betabinomial=load('beta_all.csv');
subplot(1,3,1); hold on;
P07=[kanpig(:,1),betabinomial(:,1)];
R07=[kanpig(:,2),betabinomial(:,2)];
plot(P07(:,1),R07(:,1),'.r');
plot(P07(:,2),R07(:,2),'.g');
title('All calls'); axis([0.8,1,0.6,1]); axis square; grid on; xlabel('Precision'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);

# Inside TRs
kanpig=load('kanpig_tr.csv'); 
betabinomial=load('beta_tr.csv');
subplot(1,3,2); hold on;
P07=[kanpig(:,1),betabinomial(:,1)];
R07=[kanpig(:,2),betabinomial(:,2)];
plot(P07(:,1),R07(:,1),'.r');
plot(P07(:,2),R07(:,2),'.g');
title('Inside TR'); axis([0.8,1,0.6,1]); axis square; grid on; xlabel('Precision'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);

# Outside TRs
kanpig=load('kanpig_not_tr.csv'); 
betabinomial=load('beta_not_tr.csv');
subplot(1,3,3); hold on;
P07=[kanpig(:,1),betabinomial(:,1)];
R07=[kanpig(:,2),betabinomial(:,2)];
plot(P07(:,1),R07(:,1),'.r');
plot(P07(:,2),R07(:,2),'.g');
title('Outside TR'); axis([0.8,1,0.6,1]); axis square; grid on; xlabel('Precision'); ylabel('Recall'); set(gca,'fontsize',FONT_SIZE);
