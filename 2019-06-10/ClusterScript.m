clear
% If cluster object does not exist
c = parcluster('gizmo remote R2017a');
c.AdditionalProperties.EmailAddress = 'li.lily.xie@gmail.com';
numWorkers = 24;
for i = 1:2
    cd(num2str(i));
    jobID=c.batch('main_HMSelection_mixed_par', 'Pool', numWorkers,'CurrentFolder',['/home/lxie2/shougroup/lab_users/Li/parallel/2019-06-10/' num2str(i)]);
%     jobID=c.batch('main', 'Pool', numWorkers,'CurrentFolder',['/home/lxie2/shougroup/lab_users/Li/parallel/2019-05-30/' num2str(i)]);
    cd ..
end
