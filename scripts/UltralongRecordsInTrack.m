LABELS={'TR_s','TR_e','TR_i','SG_s','SG_e','SG_i','SG10_s','SG10_e','SG10_i'};
LABELS_INS={'TR_s','SG_s','SG10_s'};

figure(1);

subplot(2,2,1); hold on;
A=load('/Users/fcunial/Downloads/ultralong_tracks/del_counts.txt');
bar(100*A([2:10],1)./A(1,1)); 
xticks(1:9); xticklabels(LABELS); ylabel('% records'); xlabel('Overlap with track'); title('DEL'); grid on; set(gca,'fontsize',14);
ylim([0 100]);

subplot(2,2,2); hold on;
A=load('/Users/fcunial/Downloads/ultralong_tracks/ins_counts.txt');
bar(100*A([2,5,8],1)./A(1,1)); 
xticks(1:3); xticklabels(LABELS_INS); ylabel('% records'); xlabel('Overlap with track'); title('INS'); grid on; set(gca,'fontsize',14);
ylim([0 100]);

subplot(2,2,3); hold on;
A=load('/Users/fcunial/Downloads/ultralong_tracks/dup_counts.txt');
bar(100*A([2:10],1)./A(1,1)); 
xticks(1:9); xticklabels(LABELS); ylabel('% records'); xlabel('Overlap with track'); title('DUP'); grid on; set(gca,'fontsize',14);
ylim([0 100]);

subplot(2,2,4); hold on;
A=load('/Users/fcunial/Downloads/ultralong_tracks/inv_counts.txt');
bar(100*A([2:10],1)./A(1,1)); 
xticks(1:9); xticklabels(LABELS); ylabel('% records'); xlabel('Overlap with track'); title('INV'); grid on; set(gca,'fontsize',14);
ylim([0 100]);
