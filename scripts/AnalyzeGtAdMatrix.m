
FONT_SIZE=16;


subplot(1,2,1); hold on;
A=load('~/Downloads/histogram_alt.txt');
plot([0:50],A(1,:)./sum(A(1,:)),'k-');
plot([0:50],A(2,:)./sum(A(2,:)),'r-');
plot([0:50],A(3,:)./sum(A(3,:)),'g-');
plot([0:50],A(4,:)./sum(A(4,:)),'b-');
axis square; grid on; title('AD ALT - After re-genotyping'); legend('all','0/0','0/1','1/1'); set(gca,'fontsize',FONT_SIZE);


subplot(1,2,2); hold on;
A=load('~/Downloads/histogram_ref.txt');
plot([0:50],A(1,:)./sum(A(1,:)),'k-');
plot([0:50],A(2,:)./sum(A(2,:)),'r-');
plot([0:50],A(3,:)./sum(A(3,:)),'g-');
plot([0:50],A(4,:)./sum(A(4,:)),'b-');
axis square; grid on; title('AD REF - After re-genotyping'); legend('all','0/0','0/1','1/1'); set(gca,'fontsize',FONT_SIZE);
