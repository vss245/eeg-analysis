%% functional connectivity tutorial - bastos 2016
% simulating data: 1 and 2 are two separate sources, 3 is the common input
% which is connected to both 1 and 2 at time delays
cfg = [];
cfg.method = 'ar';
cfg.ntrials = 200;
cfg.triallength=1;
cfg.fsample=200;
cfg.nsignal=3;
cfg.params(:,:,1)=[0.55 0 0.25;
                   0 0.55 0.25;
                   0 0 0.55]; % off-diag is the 3->1 and 3->2 influence at 1 sample
cfg.params(:,:,1)=[-0.8 0 -0.1;
                   0 -0.8 -0.1;
                   0 0 -0.8];    
cfg.noisecov = [1 0 0;0 1 0;0 0 1];
data=ft_connectivitysimulation(cfg);