SOURCE_DIR='/Users/fcunial/Downloads/BenchCohortTrios_hispanics';
FONT_SIZE=14;
CHILDREN={'1402381'};  % {'1004931','1402381','1665275','1806012','3498199'};
N_CHILDREN=length(CHILDREN);
FONT_SIZE=20;

MAX_AD=50;
MAX_LEN=10000;
MAX_DE_NOVO=0.16;

V1_ID='pileupmax_20';  #'kanpig';   % 'v1'



% % AD plot, using the old GT matrix format.
% INPUT_DIR='gt_matrix_triomerge';  % 'gt_matrix_old';
% X=[0:MAX_AD];
% for i=[1:N_CHILDREN]
%     figure(i);
%
%     subplot(2,3,1); hold on;
%     V1=load(sprintf('%s/%s/%s_%s_all_gtmatrix.txt_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
%     BB=load(sprintf('%s/%s/%s_beta_binomial_all_gtmatrix.txt_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
%     plot(X,V1,'.'); plot(X,BB,'o');
%     axis square; grid on; axis([0,MAX_AD,0,1]); xlabel('avg(AD\_REF+AD\_ALT)'); ylabel('De novo rate'); title('All calls'); set(gca,'fontsize',FONT_SIZE);
%
%     subplot(2,3,4); hold on;
%     bar(X,V1-BB);
%     axis square; grid on; axis([0,MAX_AD,-1,1]); xlabel('avg(AD\_REF+AD\_ALT)'); ylabel('V1-BB'); title('All calls'); set(gca,'fontsize',FONT_SIZE);
%
%     subplot(2,3,2); hold on;
%     V1=load(sprintf('%s/%s/%s_%s_tr_gtmatrix.txt_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
%     BB=load(sprintf('%s/%s/%s_beta_binomial_tr_gtmatrix.txt_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
%     plot(X,V1,'.'); plot(X,BB,'o');
%     axis square; grid on; axis([0,MAX_AD,0,1]); xlabel('avg(AD\_REF+AD\_ALT)'); ylabel('De novo rate'); title('Inside TR'); set(gca,'fontsize',FONT_SIZE);
%
%     subplot(2,3,5); hold on;
%     bar(X,V1-BB);
%     axis square; grid on; axis([0,MAX_AD,-1,1]); xlabel('avg(AD\_REF+AD\_ALT)'); ylabel('V1-BB'); title('Inside TR'); set(gca,'fontsize',FONT_SIZE);
%
%     subplot(2,3,3); hold on;
%     V1=load(sprintf('%s/%s/%s_%s_not_tr_gtmatrix.txt_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
%     BB=load(sprintf('%s/%s/%s_beta_binomial_not_tr_gtmatrix.txt_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
%     plot(X,V1,'.'); plot(X,BB,'o');
%     axis square; grid on; axis([0,MAX_AD,0,1]); xlabel('avg(AD\_REF+AD\_ALT)'); ylabel('De novo rate'); title('Outside TR'); set(gca,'fontsize',FONT_SIZE);
%     legend('V1','BB','location','eastoutside');
%
%     subplot(2,3,6); hold on;
%     bar(X,V1-BB);
%     axis square; grid on; axis([0,MAX_AD,-1,1]); xlabel('avg(AD\_REF+AD\_ALT)'); ylabel('V1-BB'); title('Outside TR'); set(gca,'fontsize',FONT_SIZE);
% endfor


% INS/DEL plot, using the new GT matrix format.
INPUT_DIR='.';  #'gt_matrix_triomerge';  % INPUT_DIR='gt_matrix_new';
X=[0:MAX_AD];
for i=[1:N_CHILDREN]
    figure(N_CHILDREN+i);
    
    subplot(2,3,1); hold on;
    V1=load(sprintf('%s/%s/%s_%s_all_gtmatrix.txt_ins_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_all_gtmatrix.txt_ins_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([0,MAX_AD,0,1]); xlabel('avg(AD\_REF+AD\_ALT)'); ylabel('De novo rate'); title('All INS'); set(gca,'fontsize',FONT_SIZE);
    
    subplot(2,3,2); hold on;
    V1=load(sprintf('%s/%s/%s_%s_tr_gtmatrix.txt_ins_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_tr_gtmatrix.txt_ins_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([0,MAX_AD,0,1]); xlabel('avg(AD\_REF+AD\_ALT)'); ylabel('De novo rate'); title('INS in TRs'); set(gca,'fontsize',FONT_SIZE);
    
    subplot(2,3,3); hold on;
    V1=load(sprintf('%s/%s/%s_%s_not_tr_gtmatrix.txt_ins_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_not_tr_gtmatrix.txt_ins_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([0,MAX_AD,0,1]); xlabel('avg(AD\_REF+AD\_ALT)'); ylabel('De novo rate'); title('INS outside TRs'); set(gca,'fontsize',FONT_SIZE);
    legend('V1','BB','location','eastoutside');
    
    subplot(2,3,4); hold on;
    V1=load(sprintf('%s/%s/%s_%s_all_gtmatrix.txt_del_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_all_gtmatrix.txt_del_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([0,MAX_AD,0,1]); xlabel('avg(AD\_REF+AD\_ALT)'); ylabel('De novo rate'); title('All DEL'); set(gca,'fontsize',FONT_SIZE);
    
    subplot(2,3,5); hold on;
    V1=load(sprintf('%s/%s/%s_%s_tr_gtmatrix.txt_del_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_tr_gtmatrix.txt_del_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([0,MAX_AD,0,1]); xlabel('avg(AD\_REF+AD\_ALT)'); ylabel('De novo rate'); title('DEL in TRs'); set(gca,'fontsize',FONT_SIZE);
    
    subplot(2,3,6); hold on;
    V1=load(sprintf('%s/%s/%s_%s_not_tr_gtmatrix.txt_del_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_not_tr_gtmatrix.txt_del_depth_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([0,MAX_AD,0,1]); xlabel('avg(AD\_REF+AD\_ALT)'); ylabel('De novo rate'); title('DEL outside TRs'); set(gca,'fontsize',FONT_SIZE);
    % legend('V1','BB','location','eastoutside');
endfor


% SVLEN plot, using the new GT matrix format.
INPUT_DIR='.';  #'gt_matrix_triomerge';  % INPUT_DIR='gt_matrix_new';
X=[20,50,100,200,300,400,500,600,700,800,1000,2000,3000,4000,5000,6000,7000,8000,9000,MAX_LEN];
for i=[1:N_CHILDREN]
    figure(2*N_CHILDREN+i);
    
    subplot(2,3,1); hold on;
    V1=load(sprintf('%s/%s/%s_%s_all_gtmatrix.txt_ins_len_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_all_gtmatrix.txt_ins_len_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([0,1000,0,0.6]); xlabel('SVLEN'); ylabel('De novo rate'); title('All INS'); set(gca,'fontsize',FONT_SIZE);
    
    subplot(2,3,2); hold on;
    V1=load(sprintf('%s/%s/%s_%s_tr_gtmatrix.txt_ins_len_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_tr_gtmatrix.txt_ins_len_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([0,1000,0,0.6]); xlabel('SVLEN'); ylabel('De novo rate'); title('INS in TRs'); set(gca,'fontsize',FONT_SIZE);
    
    subplot(2,3,3); hold on;
    V1=load(sprintf('%s/%s/%s_%s_not_tr_gtmatrix.txt_ins_len_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_not_tr_gtmatrix.txt_ins_len_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([0,1000,0,0.6]); xlabel('SVLEN'); ylabel('De novo rate'); title('INS outside TRs'); set(gca,'fontsize',FONT_SIZE);
    legend('V1','BB','location','eastoutside');
    
    subplot(2,3,4); hold on;
    V1=load(sprintf('%s/%s/%s_%s_all_gtmatrix.txt_del_len_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_all_gtmatrix.txt_del_len_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([0,1000,0,0.6]); xlabel('SVLEN'); ylabel('De novo rate'); title('All DEL'); set(gca,'fontsize',FONT_SIZE);
    
    subplot(2,3,5); hold on;
    V1=load(sprintf('%s/%s/%s_%s_tr_gtmatrix.txt_del_len_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_tr_gtmatrix.txt_del_len_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([0,1000,0,0.6]); xlabel('SVLEN'); ylabel('De novo rate'); title('DEL in TRs'); set(gca,'fontsize',FONT_SIZE);
    
    subplot(2,3,6); hold on;
    V1=load(sprintf('%s/%s/%s_%s_not_tr_gtmatrix.txt_del_len_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_not_tr_gtmatrix.txt_del_len_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([0,1000,0,0.6]); xlabel('SVLEN'); ylabel('De novo rate'); title('DEL outside TRs'); set(gca,'fontsize',FONT_SIZE);
    % legend('V1','BB','location','eastoutside');
endfor


% NUMNEIGH plot, using the new GT matrix format.
INPUT_DIR='.';  #'gt_matrix_triomerge';  % INPUT_DIR='gt_matrix_numneigh';
X=[0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000];
for i=[1:N_CHILDREN]
    figure(3*N_CHILDREN+i);
    
    subplot(2,3,1); hold on;
    V1=load(sprintf('%s/%s/%s_%s_all_gtmatrix.txt_numneigh_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_all_gtmatrix.txt_numneigh_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([-1,11,0,0.6]); xlabel('NumNeighbors'); ylabel('De novo rate'); title('All calls'); set(gca,'fontsize',FONT_SIZE);
    
    subplot(2,3,4); hold on;
    V1=load(sprintf('%s/%s/%s_%s_all_gtmatrix.txt_numneigh.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_all_gtmatrix.txt_numneigh.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; xlabel('NumNeighbors'); ylabel('Number of calls'); title('All calls'); set(gca,'fontsize',FONT_SIZE); axis([-1,11,0,12000]); 
    
    subplot(2,3,2); hold on;
    V1=load(sprintf('%s/%s/%s_%s_tr_gtmatrix.txt_numneigh_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_tr_gtmatrix.txt_numneigh_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([-1,11,0,0.6]); xlabel('NumNeighbors'); ylabel('De novo rate'); title('Inside TRs'); set(gca,'fontsize',FONT_SIZE);
    
    subplot(2,3,5); hold on;
    V1=load(sprintf('%s/%s/%s_%s_tr_gtmatrix.txt_numneigh.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_tr_gtmatrix.txt_numneigh.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; xlabel('NumNeighbors'); ylabel('Number of calls'); title('Inside TRs'); set(gca,'fontsize',FONT_SIZE); axis([-1,11,0,12000]); 
    
    subplot(2,3,3); hold on;
    V1=load(sprintf('%s/%s/%s_%s_not_tr_gtmatrix.txt_numneigh_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_not_tr_gtmatrix.txt_numneigh_to_denovo_rate.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; axis([-1,11,0,0.6]); xlabel('NumNeighbors'); ylabel('De novo rate'); title('Outside TRs'); set(gca,'fontsize',FONT_SIZE);
    legend('V1','BB','location','eastoutside');
    
    subplot(2,3,6); hold on;
    V1=load(sprintf('%s/%s/%s_%s_not_tr_gtmatrix.txt_numneigh.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i},V1_ID));
    % BB=load(sprintf('%s/%s/%s_beta_binomial_not_tr_gtmatrix.txt_numneigh.txt',SOURCE_DIR,INPUT_DIR,CHILDREN{i}));
    plot(X,V1,'.'); 
    % plot(X,BB,'o');
    axis square; grid on; xlabel('NumNeighbors'); ylabel('Number of calls'); title('Outside TRs'); set(gca,'fontsize',FONT_SIZE); axis([-1,11,0,12000]);
    % legend('V1','BB','location','eastoutside');
endfor
