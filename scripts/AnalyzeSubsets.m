SOURCE_DIR='~/Downloads';
MAX_AD=30;
FONT_SIZE=10;





A=load(sprintf('%s/ad_all_cleaned.txt',SOURCE_DIR));
[nrows,ncolumns]=size(A);

subplot(3,4,1);
plot(A(:,1),A(:,2),'.');
xlabel('AD REF'); ylabel('AD ALT'); title('Sample 1328181 - All calls'); set(gca,'fontsize',FONT_SIZE); axis square; grid on;
B=A(:,1)+A(:,2)>MAX_AD;
printf('Fraction of records with AD_REF+AD_ALT>%d: %f %%\n', MAX_AD, 100 * nnz(B) / nrows);

subplot(3,4,2);
plot(A(:,1),A(:,2),'.');
xlabel('AD REF'); ylabel('AD ALT'); title(''); set(gca,'fontsize',FONT_SIZE); axis square; grid on;
axis([0,150,0,150]);

subplot(3,4,3);
matrix=zeros(MAX_AD+1,MAX_AD+1);
for i=[1:nrows]
    x=min(A(i,1),MAX_AD);
    y=min(A(i,2),MAX_AD);
    matrix(x+1,y+1)=matrix(x+1,y+1)+1;
endfor
imagesc(matrix); colormap('hot'); xlabel('AD REF'); ylabel('AD ALT'); title(''); axis square; set(gca,'fontsize',FONT_SIZE); 

subplot(3,4,4);
B=A(:,1)+A(:,2);
[Y,X]=hist(B,1000); plot(X,log10(Y),'.'); axis square; grid on; xlabel('AD REF + AD ALT'); ylabel('n. records (log10)'); set(gca,'fontsize',FONT_SIZE); 





A=load(sprintf('%s/ad_present_cleaned.txt',SOURCE_DIR));
[nrows,ncolumns]=size(A);

subplot(3,4,5);
plot(A(:,1),A(:,2),'.');
xlabel('AD REF'); ylabel('AD ALT'); title('Sample 1328181 - Present calls'); set(gca,'fontsize',FONT_SIZE); axis square; grid on;
B=A(:,1)+A(:,2)>MAX_AD;
printf('Fraction of records with AD_REF+AD_ALT>%d: %f %%\n', MAX_AD, 100 * nnz(B) / nrows);

subplot(3,4,6);
plot(A(:,1),A(:,2),'.');
xlabel('AD REF'); ylabel('AD ALT'); title(''); set(gca,'fontsize',FONT_SIZE); axis square; grid on;
axis([0,150,0,150]);

subplot(3,4,7);
matrix=zeros(MAX_AD+1,MAX_AD+1);
for i=[1:nrows]
    x=min(A(i,1),MAX_AD);
    y=min(A(i,2),MAX_AD);
    matrix(x+1,y+1)=matrix(x+1,y+1)+1;
endfor
imagesc(matrix); colormap('hot'); xlabel('AD REF'); ylabel('AD ALT'); title(''); axis square; set(gca,'fontsize',FONT_SIZE); 

subplot(3,4,8);
B=A(:,1)+A(:,2);
[Y,X]=hist(B,1000); plot(X,log10(Y),'.'); axis square; grid on; xlabel('AD REF + AD ALT'); ylabel('n. records (log10)'); set(gca,'fontsize',FONT_SIZE); 




A=load(sprintf('%s/ad_present_and_50_cleaned.txt',SOURCE_DIR));
[nrows,ncolumns]=size(A);

subplot(3,4,9);
plot(A(:,1),A(:,2),'.');
xlabel('AD REF'); ylabel('AD ALT'); title('Sample 1328181 - Present calls >=50bp'); set(gca,'fontsize',FONT_SIZE); axis square; grid on;
B=A(:,1)+A(:,2)>MAX_AD;
printf('Fraction of records with AD_REF+AD_ALT>%d: %f %%\n', MAX_AD, 100 * nnz(B) / nrows);

subplot(3,4,10);
plot(A(:,1),A(:,2),'.');
xlabel('AD REF'); ylabel('AD ALT'); title(''); set(gca,'fontsize',FONT_SIZE); axis square; grid on;
axis([0,150,0,150]);

subplot(3,4,11);
matrix=zeros(MAX_AD+1,MAX_AD+1);
for i=[1:nrows]
    x=min(A(i,1),MAX_AD);
    y=min(A(i,2),MAX_AD);
    matrix(x+1,y+1)=matrix(x+1,y+1)+1;
endfor
imagesc(matrix); colormap('hot'); xlabel('AD REF'); ylabel('AD ALT'); title(''); axis square; set(gca,'fontsize',FONT_SIZE); 

subplot(3,4,12);
B=A(:,1)+A(:,2);
[Y,X]=hist(B,1000); plot(X,log10(Y),'.'); axis square; grid on; xlabel('AD REF + AD ALT'); ylabel('n. records (log10)'); set(gca,'fontsize',FONT_SIZE); 
