X=[0:100];
N_SAMPLES=[128];

for i=N_SAMPLES
	figure(i);
	A=load('~/Downloads/HG00097_kanpig_annotated_histogram.txt');
	
	subplot(1,2,1); hold on;
	plot(X,A(:,1)./sum(A(:,1)),'.r'); plot(X,A(:,4)./sum(A(:,4)),'.b');
	axis square; grid on; xlabel('KS'); ylabel('Fraction of all calls of each category'); title(sprintf('>=%d samples / HG00097 / Calls with just one KS value',i)); legend('TP','FP','location','northwest'); set(gca,'fontsize',16);
	%axis([80,100,0,1400]);

	subplot(1,2,2); hold on;
	plot(X,A(:,2)./sum(A(:,2)),'.r'); plot(X,A(:,3)./sum(A(:,3)),'.m'); plot(X,A(:,5)./sum(A(:,5)),'.b'); plot(X,A(:,6)./sum(A(:,6)),'.c');
	axis square; grid on; xlabel('KS1 and KS2'); ylabel('Fraction of all calls of each category'); title(sprintf('>=%d samples / HG00097 / Calls with two KS values',i)); legend('TP KS1','TP KS2','FP KS1','FP KS2','location','northwest'); set(gca,'fontsize',16);
	%axis([60,100,0,600]);
end
