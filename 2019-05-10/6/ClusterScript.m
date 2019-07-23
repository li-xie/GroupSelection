clear
% If cluster object does not exist
c = parcluster;
% c.AdditionalProperties.DebugMessagesTurnedOn=true;
c.AdditionalProperties.AdditionalSubmitArgs= '--exclude=gizmof287,gizmof288';
c.AdditionalProperties.EmailAddress = 'aeyuan@uw.edu';
numWorkers = 49;
jobID=c.batch('main', 'Pool', numWorkers,'CurrentFolder','/home/ayuan/shougroup/lab_users/Alex_Yuan/projects/community_selection/code_group_selection/');
