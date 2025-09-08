SOURCE_DIR='/Users/fcunial/Downloads/gt_priors/truvari_out';
FONT_SIZE=14;

AD_MAX=50;

figure(1);

subplot(3,3,1);
A=load(sprintf('%s/distributions_00.csv',SOURCE_DIR));
B=A./sum(A,2);
imagesc(B); colormap('hot'); axis([0,20,0,20]); axis square; xlabel('AD\_ALT'); ylabel('Depth'); title('0/0 - AD\_ALT distribution at each depth'); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,4); hold on;
plot([0:AD_MAX]./5,B(6,:)'); 
plot([0:AD_MAX]./10,B(11,:)');
axis([0,1,0,0.6]); axis square; grid on; xlabel('Fraction of depth'); ylabel('p(AD\_ALT)'); title('Zoom in on two depths'); legend({'Depth 5','Depth 10'}); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,7); hold on;
S=zeros(AD_MAX+1,1);
for i=[0:AD_MAX]
    for j=[0:AD_MAX]
        S(i+1)=S(i+1)+A(i+1,j+1)*(j/i);
    endfor
    S(i+1)=S(i+1)/sum(A(i+1,:));
endfor
plot([0:AD_MAX],S,'.'); axis([0,30,0,1]); axis square; grid on; xlabel('Depth'); title('S(0/0)'); set(gca,'fontsize',FONT_SIZE);




subplot(3,3,2);
A=load(sprintf('%s/distributions_01.csv',SOURCE_DIR));
B=A./sum(A,2);
imagesc(B); colormap('hot'); axis([0,20,0,20]); axis square; xlabel('AD\_ALT'); ylabel('Depth'); title('0/1 - AD\_ALT distribution at each depth'); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,5); hold on;
plot([0:AD_MAX]./5,B(6,:)'); 
plot([0:AD_MAX]./10,B(11,:)');
axis([0,1,0,0.6]); axis square; grid on; xlabel('Fraction of depth'); ylabel('p(AD\_ALT)'); title('Zoom in on two depths'); legend({'Depth 5','Depth 10'}); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,8); hold on;
S=zeros(AD_MAX+1,1);
for i=[0:AD_MAX]
    for j=[0:AD_MAX]
        S(i+1)=S(i+1)+A(i+1,j+1)*(j/i);
    endfor
    S(i+1)=S(i+1)/sum(A(i+1,:));
endfor
plot([0:AD_MAX],S,'.'); axis([0,30,0,1]); axis square; grid on; xlabel('Depth'); title('S(0/1)'); set(gca,'fontsize',FONT_SIZE);




subplot(3,3,3);
A=load(sprintf('%s/distributions_11.csv',SOURCE_DIR));
B=A./sum(A,2);
imagesc(B); colormap('hot'); axis([0,20,0,20]); axis square; xlabel('AD\_ALT'); ylabel('Depth'); title('1/1 - AD\_ALT distribution at each depth'); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,6); hold on;
plot([0:AD_MAX]./5,B(6,:)'); 
plot([0:AD_MAX]./10,B(11,:)');
axis([0,1,0,1]); axis square; grid on; xlabel('Fraction of depth'); ylabel('p(AD\_ALT)'); title('Zoom in on two depths'); legend({'Depth 5','Depth 10'},'location','northwest'); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,9); hold on;
S=zeros(AD_MAX+1,1);
for i=[0:AD_MAX]
    for j=[0:AD_MAX]
        S(i+1)=S(i+1)+A(i+1,j+1)*(j/i);
    endfor
    S(i+1)=S(i+1)/sum(A(i+1,:));
endfor
plot([0:AD_MAX],S,'.'); axis([0,30,0,1]); axis square; grid on; xlabel('Depth'); title('S(1/1)'); set(gca,'fontsize',FONT_SIZE);
