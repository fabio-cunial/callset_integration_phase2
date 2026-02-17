
TITLES={'Records in chr','Records in chr without FILTER','Records in BED','Records in BED without FILTER'};

figure(1);
A=load('raw_data/matrix.csv');
for i=[1:4]
    subplot(4,1,i);
    [yy,xx]=hist(A(:,i));
    bar(xx,yy);
    axis square; ylabel('n. samples'); xlabel('n. calls'); title(TITLES{i});
    %grid on; set(gca,'fontsize',12);
endfor
