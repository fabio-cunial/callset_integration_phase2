
subplot(1,4,1); A=load('matrix_del.csv'); imagesc(log10(A)); colormap hot; ylabel('N\_DISCOVERY\_SAMPLES'); xlabel('SVLEN (log10)'); title('Ultralong DEL'); set(gca,'fontsize',14); axis([4.5,9.5,1,700]);
subplot(1,4,2); A=load('matrix_ins.csv'); imagesc(log10(A)); colormap hot; ylabel('N\_DISCOVERY\_SAMPLES'); xlabel('SVLEN (log10)'); title('Ultralong INS'); set(gca,'fontsize',14); axis([4.5,9.5,1,700]);
subplot(1,4,3); A=load('matrix_dup.csv'); imagesc(log10(A)); colormap hot; ylabel('N\_DISCOVERY\_SAMPLES'); xlabel('SVLEN (log10)'); title('Ultralong DUP'); set(gca,'fontsize',14); axis([4.5,9.5,1,700]);
subplot(1,4,4); A=load('matrix_inv.csv'); imagesc(log10(A)); colormap hot; ylabel('N\_DISCOVERY\_SAMPLES'); xlabel('SVLEN (log10)'); title('Ultralong INV'); set(gca,'fontsize',14); axis([4.5,9.5,1,700]);
