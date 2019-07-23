% clear
% close all
% top-dog
% comm_f = [3 11 12];
% group_f = [4 13 14];
% % top-10
% comm_f = [5 7 8];
% group_f = [6 9 10];

% figure(2)
% % hold on
% % for i = 1:3
% %     load([num2str(comm_f(i)) '/fp_data']);
% %     plot(log10(mean(varfp0,2)),'k')
% % end
% % 
for i = 1:3
    load([num2str((i+3)) '/fp_data']);
    plot(log10(mean(varfp0,2)),'r')   
end
hold off
xlabel('Cycle')
ylabel('log10 selected var(0)')
axis([1 1e3 -8 -2])

% figure(1)
% hold on
% for i = 1:3
%     load([num2str(comm_f(i)) '/fp_data']);
%     plot(mean(fp0,2),'k')
% end

% for i = 1:3
%     load([num2str(i+3) '/fp_data']);
%     plot(mean(fp0,2),'r')   
% end
% ylim([0.11 0.14])
% hold off
% xlabel('Cycle')
% ylabel('selected fp(0)')