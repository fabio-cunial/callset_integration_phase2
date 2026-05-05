A=load('~/Downloads/del_cleaned.csv');
subplot(1,2,1); hold on; plot(A(:,1)./sum(A(:,1)),'.b');
subplot(1,2,2); hold on; plot(A(:,2)./sum(A(:,2)),'.b');

A=load('~/Downloads/ins_cleaned.csv');
subplot(1,2,1); hold on; plot(A(:,1)./sum(A(:,1)),'.r');
subplot(1,2,2); hold on; plot(A(:,2)./sum(A(:,2)),'.r');

A=load('~/Downloads/dup_cleaned.csv');
subplot(1,2,1); hold on; plot(A(:,1)./sum(A(:,1)),'.m');
subplot(1,2,2); hold on; plot(A(:,2)./sum(A(:,2)),'.m');

A=load('~/Downloads/inv_cleaned.csv');
subplot(1,2,1); hold on; plot(A(:,1)./sum(A(:,1)),'.g');
subplot(1,2,2); hold on; plot(A(:,2)./sum(A(:,2)),'.g');

subplot(1,2,1); axis square; xlabel('chr'); ylabel('Fraction of calls'); title('Ultralong calls, total'); legend('DEL','INS','DUP','INV'); grid on; set(gca,'fontsize',18);
subplot(1,2,2); axis square; xlabel('chr'); ylabel('Fraction of calls'); title('Ultralong calls, true'); legend('DEL','INS','DUP','INV'); grid on; set(gca,'fontsize',18);
