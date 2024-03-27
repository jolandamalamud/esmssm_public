addpath(genpath('/home/skgtjmw/Scratch/code/cluster/'))

% Get the Task ID of this job
job_id_string = getenv('SGE_TASK_ID');
con_in = str2double(job_id_string);

option.filepath = '/home/skgtjmw/Scratch/ssm';
option.foldername = 'twindataset';%'simulations'; %'twindataset'; 'leemput_patients'; 'leemput_controls'; 'esm_merge'
option.goalpath = ['/home/skgtjmw/Scratch/ssm/data_fit/' option.foldername];
option.input = true;
option.bias = true;
option.maxms = 10;
option.r = 1000;
option.k = 60;
option.zero_centered = false; % if data is prepocessed to be centered around 0

option.moodfiles = dir(sprintf('%s/data_prep_cluster/%s/datafit_*.mat', ...
            option.filepath, option.foldername));
                
Nsj = length(option.moodfiles);

subj_list = [1:Nsj];

id = subj_list(con_in);

% This script runs the jobs
disp(['starting ssm_job.m for id: ' mat2str(id)])

ssm_emfit_fast_input(id, option);

% If we have a job ID, then exit so matlab does not sit at the command
% prompt forever.
if con_in ~= -1
    exit
end
