clear
% If cluster object does not exist
c = parcluster('gizmo remote R2017a');
c.AdditionalProperties.EmailAddress = 'li.lily.xie@gmail.com';
numWorkers = 24;
for i=5:6
    cd(num2str(i));
    jobID=c.batch('main', 'Pool', numWorkers,'CurrentFolder',['/home/lxie2/shougroup/lab_users/Li/parallel/2019-05-10/' num2str(i)]);
    cd ..
end
