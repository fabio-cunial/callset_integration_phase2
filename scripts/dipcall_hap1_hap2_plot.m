FONT_SIZE=18;

figure(1);

subplot(1,2,1); hold on;
A=load('coverage_hap1.csv');
B=load('coverage_hap2.csv');
plot(B(:,1),B(:,2),'or');
plot(A(:,1),A(:,2),'.b');
xlabel('chrX'); ylabel('chrY'); title('Percentage of covered bases'); axis square; axis([0,100,0,100]); grid on; set(gca,'fontsize',FONT_SIZE);

subplot(1,2,2); hold on;
A=load('meandepth_hap1.csv');
B=load('meandepth_hap2.csv');
plot(B(:,1),B(:,2),'or');
plot(A(:,1),A(:,2),'.b');
xlabel('chrX'); ylabel('chrY'); title('Mean depth of coverage'); axis square; axis([0,1,0,1]); grid on; set(gca,'fontsize',FONT_SIZE);
