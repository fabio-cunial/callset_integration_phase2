SOURCE_DIR='/Users/fcunial/Downloads/gt_priors/truvari_stringent';
FONT_SIZE=14;
SAMPLES={"HG002","HG00438","HG005","HG00621","HG00673","HG00733","HG00741","HG01071","HG01106","HG01109","HG01123","HG01175","HG01243","HG01258","HG01358","HG01361","HG01891","HG01928","HG01952","HG01978","HG02055","HG02080","HG02109","HG02145","HG02148","HG02257","HG02486","HG02559","HG02572","HG02622","HG02630","HG02717","HG02723","HG02818","HG02886","HG03098","HG03453","HG03486","HG03492","HG03516","HG03540","HG03579","NA18906","NA19240","NA20129"};
printf("Considering only %d samples\n",length(SAMPLES));

AD_MAX=50;







figure(1);



subplot(3,3,1);
A=zeros(AD_MAX+1,AD_MAX+1);
for i=[1:length(SAMPLES)]
    B=load(sprintf('%s/%s_distributions_00.csv',SOURCE_DIR,SAMPLES{i}));
    A=A+B;
endfor
B=A./sum(A,2);
imagesc(B); colormap('hot'); axis([0,20,0,20]); axis square; xlabel('AD\_ALT'); ylabel('Depth'); title('0/0 - AD\_ALT distribution at each depth'); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,4); hold on;
plot([0:AD_MAX]./5,B(6,:)'); 
plot([0:AD_MAX]./10,B(11,:)');
axis([0,1,0,0.6]); axis square; grid on; xlabel('Fraction of depth'); ylabel('p(AD\_ALT)'); title('Zoom in on two depths'); legend({'Depth 5','Depth 10'}); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,7); hold on;
S=zeros(AD_MAX+1,1);
for i=[0:AD_MAX]
    for j=[0:i]
        S(i+1)=S(i+1)+A(i+1,j+1)*(j/i);
    endfor
    S(i+1)=S(i+1)/sum(A(i+1,:));
endfor
plot([0:AD_MAX],S,'.'); axis([0,30,0,1]); axis square; grid on; xlabel('Depth'); title('S(0/0)'); set(gca,'fontsize',FONT_SIZE);




subplot(3,3,2);
A=zeros(AD_MAX+1,AD_MAX+1);
for i=[1:length(SAMPLES)]
    B=load(sprintf('%s/%s_distributions_01.csv',SOURCE_DIR,SAMPLES{i}));
    A=A+B;
endfor
B=A./sum(A,2);
imagesc(B); colormap('hot'); axis([0,20,0,20]); axis square; xlabel('AD\_ALT'); ylabel('Depth'); title('0/1 - AD\_ALT distribution at each depth'); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,5); hold on;
plot([0:AD_MAX]./5,B(6,:)'); 
plot([0:AD_MAX]./10,B(11,:)');
axis([0,1,0,0.6]); axis square; grid on; xlabel('Fraction of depth'); ylabel('p(AD\_ALT)'); title('Zoom in on two depths'); legend({'Depth 5','Depth 10'}); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,8); hold on;
S=zeros(AD_MAX+1,1);
for i=[0:AD_MAX]
    for j=[0:i]
        S(i+1)=S(i+1)+A(i+1,j+1)*(j/i);
    endfor
    S(i+1)=S(i+1)/sum(A(i+1,:));
endfor
plot([0:AD_MAX],S,'.'); axis([0,30,0,1]); axis square; grid on; xlabel('Depth'); title('S(0/1)'); set(gca,'fontsize',FONT_SIZE);




subplot(3,3,3);
A=zeros(AD_MAX+1,AD_MAX+1);
for i=[1:length(SAMPLES)]
    B=load(sprintf('%s/%s_distributions_11.csv',SOURCE_DIR,SAMPLES{i}));
    A=A+B;
endfor
B=A./sum(A,2);
imagesc(B); colormap('hot'); axis([0,20,0,20]); axis square; xlabel('AD\_ALT'); ylabel('Depth'); title('1/1 - AD\_ALT distribution at each depth'); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,6); hold on;
plot([0:AD_MAX]./5,B(6,:)'); 
plot([0:AD_MAX]./10,B(11,:)');
axis([0,1,0,1]); axis square; grid on; xlabel('Fraction of depth'); ylabel('p(AD\_ALT)'); title('Zoom in on two depths'); legend({'Depth 5','Depth 10'},'location','northwest'); set(gca,'fontsize',FONT_SIZE);

subplot(3,3,9); hold on;
S=zeros(AD_MAX+1,1);
for i=[0:AD_MAX]
    for j=[0:i]
        S(i+1)=S(i+1)+A(i+1,j+1)*(j/i);
    endfor
    S(i+1)=S(i+1)/sum(A(i+1,:));
endfor
plot([0:AD_MAX],S,'.'); axis([0,30,0,1]); axis square; grid on; xlabel('Depth'); title('S(1/1)'); set(gca,'fontsize',FONT_SIZE);
