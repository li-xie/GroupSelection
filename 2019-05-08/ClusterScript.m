clear
% If cluster object does not exist
c = parcluster('gizmo remote R2017a');
c.AdditionalProperties.EmailAddress = 'li.lily.xie@gmail.com';
numWorkers = 24;
for i=13:14
    cd(num2str(i));
%     jobID=c.batch('main_HMZeroNull_par', 'Pool', numWorkers,'CurrentFolder',['/home/lxie2/shougroup/lab_users/Li/parallel/2019-05-08/' num2str(i)]);
    jobID=c.batch('main', 'Pool', numWorkers,'CurrentFolder',['/home/lxie2/shougroup/lab_users/Li/parallel/2019-05-08/' num2str(i)]);
    cd ..
end
