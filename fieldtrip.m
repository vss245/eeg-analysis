%% Fieldtrip preprocessing pipeline
cfg = [];
cfg.dataset     = 'subj2.vhdr';
data_eeg        = ft_preprocessing(cfg);