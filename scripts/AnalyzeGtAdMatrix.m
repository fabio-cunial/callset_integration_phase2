INPUT_DIR='~/Downloads/AnalyzeGtAdMatrix/new';
FONT_SIZE=16;
MAX_AD=100;
QUANTUM=0.01;




% Counts
figure(1);
subplot(2,2,1); hold on;
A=load(sprintf('%s/before_regenotyping/histogram_alt.txt',INPUT_DIR));
plot([0:MAX_AD],A(2,:)./sum(A(2,:)),'r-');
plot([0:MAX_AD],A(3,:)./sum(A(3,:)),'g-');
plot([0:MAX_AD],A(4,:)./sum(A(4,:)),'b-');
axis([0,30,0,0.15]); axis square; grid on; title('AD ALT - Before re-genotyping'); ylabel('% of calls'); legend('0/0','0/1','1/1'); set(gca,'fontsize',FONT_SIZE);

subplot(2,2,2); hold on;
A=load(sprintf('%s/before_regenotyping/histogram_ref.txt',INPUT_DIR));
plot([0:MAX_AD],A(2,:)./sum(A(2,:)),'r-');
plot([0:MAX_AD],A(3,:)./sum(A(3,:)),'g-');
plot([0:MAX_AD],A(4,:)./sum(A(4,:)),'b-');
axis([0,30,0,0.15]); axis square; grid on; title('AD REF - Before re-genotyping'); ylabel('% of calls'); legend('0/0','0/1','1/1'); set(gca,'fontsize',FONT_SIZE);

subplot(2,2,3); hold on;
A=load(sprintf('%s/after_regenotyping/histogram_alt.txt',INPUT_DIR));
plot([0:MAX_AD],A(2,:)./sum(A(2,:)),'r-');
plot([0:MAX_AD],A(3,:)./sum(A(3,:)),'g-');
plot([0:MAX_AD],A(4,:)./sum(A(4,:)),'b-');
axis([0,30,0,0.15]); axis square; grid on; title('AD ALT - After re-genotyping'); ylabel('% of calls'); legend('0/0','0/1','1/1'); set(gca,'fontsize',FONT_SIZE);

subplot(2,2,4); hold on;
A=load(sprintf('%s/after_regenotyping/histogram_ref.txt',INPUT_DIR));
plot([0:MAX_AD],A(2,:)./sum(A(2,:)),'r-');
plot([0:MAX_AD],A(3,:)./sum(A(3,:)),'g-');
plot([0:MAX_AD],A(4,:)./sum(A(4,:)),'b-');
axis([0,30,0,0.15]); axis square; grid on; title('AD REF - After re-genotyping'); ylabel('% of calls'); legend('0/0','0/1','1/1'); set(gca,'fontsize',FONT_SIZE);




% Fractions
figure(2);
subplot(2,2,1); hold on;
A=load(sprintf('%s/before_regenotyping/histogram_alt_fractional.txt',INPUT_DIR));
plot([0:QUANTUM:1],A(2,:)./sum(A(2,:)),'r-');
plot([0:QUANTUM:1],A(3,:)./sum(A(3,:)),'g-');
plot([0:QUANTUM:1],A(4,:)./sum(A(4,:)),'b-');
axis([0,1,0,0.1]); axis square; grid on; title('AD ALT - Before re-genotyping'); ylabel('% of calls'); legend('0/0','0/1','1/1'); set(gca,'fontsize',FONT_SIZE);

subplot(2,2,2); hold on;
A=load(sprintf('%s/before_regenotyping/histogram_ref_fractional.txt',INPUT_DIR));
plot([0:QUANTUM:1],A(2,:)./sum(A(2,:)),'r-');
plot([0:QUANTUM:1],A(3,:)./sum(A(3,:)),'g-');
plot([0:QUANTUM:1],A(4,:)./sum(A(4,:)),'b-');
axis([0,1,0,0.1]); axis square; grid on; title('AD REF - Before re-genotyping'); ylabel('% of calls'); legend('0/0','0/1','1/1'); set(gca,'fontsize',FONT_SIZE);

subplot(2,2,3); hold on;
A=load(sprintf('%s/after_regenotyping/histogram_alt_fractional.txt',INPUT_DIR));
plot([0:QUANTUM:1],A(2,:)./sum(A(2,:)),'r-');
plot([0:QUANTUM:1],A(3,:)./sum(A(3,:)),'g-');
plot([0:QUANTUM:1],A(4,:)./sum(A(4,:)),'b-');
axis([0,1,0,0.1]); axis square; grid on; title('AD ALT - After re-genotyping'); ylabel('% of calls'); legend('0/0','0/1','1/1'); set(gca,'fontsize',FONT_SIZE);

subplot(2,2,4); hold on;
A=load(sprintf('%s/after_regenotyping/histogram_ref_fractional.txt',INPUT_DIR));
plot([0:QUANTUM:1],A(2,:)./sum(A(2,:)),'r-');
plot([0:QUANTUM:1],A(3,:)./sum(A(3,:)),'g-');
plot([0:QUANTUM:1],A(4,:)./sum(A(4,:)),'b-');
axis([0,1,0,0.1]); axis square; grid on; title('AD REF - After re-genotyping'); ylabel('% of calls'); legend('0/0','0/1','1/1'); set(gca,'fontsize',FONT_SIZE);
