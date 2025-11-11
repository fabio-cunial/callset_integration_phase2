X=[0:100];
N_SAMPLES=[1,2,32,1024,2048];

for i=N_SAMPLES
	figure(i);
	A=load(sprintf('~/Downloads/1000920_%d_samples_histogram.txt',i));
	
	subplot(1,2,1); hold on;
	plot(X,A(:,1),'.r'); plot(X,A(:,4),'.b');
	axis square; grid on; xlabel('KS'); title(sprintf('>=%d samples / 1000920 / Calls with just one KS value',i)); legend('present','absent','location','northwest'); set(gca,'fontsize',16);
	%axis([80,100,0,1400]);

	subplot(1,2,2); hold on;
	plot(X,A(:,2),'.r'); plot(X,A(:,3),'.m'); plot(X,A(:,5),'.b'); plot(X,A(:,6),'.c');
	axis square; grid on; xlabel('KS1 and KS2'); title(sprintf('>=%d samples / 1000920 / Calls with two KS values',i)); legend('present KS1','present KS2','absent KS1','absent KS2','location','northwest'); set(gca,'fontsize',16);
	%axis([60,100,0,600]);
end
