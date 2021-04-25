% Basic preprocessing pipeline in FieldTrip (based on CamCAN data)
%% Load the data for one subject
% Paths are set elsewhere
function preprocess(nsub)
  cfg = [];
  subj = ['sub-CC' num2str(sub_info(nsub))];
  save_prep = [preprocessed_dir];
  if (isfolder(data_dir))
      filename = [[meg_dir subj] '/ses-rest/meg/' subj '_ses-rest_task-rest_proc-sss.fif'];
      cfg.dataset = filename; %for fif data
      % Settings
      cfg.continuous='yes';
      cfg.channel='all';
      cfg.demean='yes';
      cfg.bpfilter = 'yes'; % Bandpass filter
      cfg.bpfreq = [1 45];
      cfg.bsfilter = 'yes'; % Bandstop to remove line noise (~50 Hz)
      cfg.bsfreq = [48 52];
      data=ft_preprocessing(cfg);
      % Get layout for plots
      grad = ft_read_sens(filename, 'senstype', 'meg');
      cfg = [];
      cfg.grad = grad;
      cfg.projection = 'polar';
      cfg.channel = 'MEG*1'; % these are the magnetometers
      layout.mag = ft_prepare_layout(cfg);
      cfg.channel = 'MEG*2'; % planar gradiometers
      layout.grad1 = ft_prepare_layout(cfg);
      cfg.channel = 'MEG*3'; % planar gradiometers
      layout.grad2 = ft_prepare_layout(cfg);
      % Find indices for channels (and practice writing regexes)
      eogs = find(~cellfun(@isempty,regexp(data.label,'^EOG')));
      ecg = find(~cellfun(@isempty,regexp(data.label,'^ECG')));
      chantypes(:,1) = find(~cellfun(@isempty,regexp(data.label,'MEG\d*1$')));
      chantypes(:,2) = find(~cellfun(@isempty,regexp(data.label,'MEG\d*2$')));
      chantypes(:,3) = find(~cellfun(@isempty,regexp(data.label,'MEG\d*3$')));
      % Downsample to 100 Hz
      cfg = [];
      cfg.resamplefs = 100;
      cfg.method = 'downsample';
      data=ft_resampledata(cfg,data);
      % Regress out EOG and ECG data
      eeg_data = data.trial{1,1}';
      eog_data = data.trial{1,1}(eogs,:)';
      ecg_data = data.trial{1,1}(ecg,:)';
      eeg_data = eeg_data - eog_data*(eog_data\eeg_data);
      eeg_data = eeg_data - ecg_data*(ecg_data\eeg_data);
      data.trial{1,1} = eeg_data';
      % Reject outlying channels (3 scaled median deviations from the median)
      % and interpolate them using a weighted average of neighbors
      % (implemented using a helper function, can also be done in FT but manually)
      chanlist = 1:length(data.label);
      all_bad_ch = [];
      for chtype = 1:size(chantypes,2)
          [interdata, bc] = detect_bad_channels(save_report, chantypes(:,chtype), data, layout.mag);
          all_bad_ch = [all_bad_ch;bc];
      end
      % Cut the data into 2 second segments
      cfg = [];
      cfg.length = 2;
      cfg.overlap = 0;
      data_seg = ft_redefinetrial(cfg,interdata);
      % Repeat outlier rejection, but for trials
      all_bad_tri = [];
      for chtype = 1:size(chantypes,2)
          bt = detect_bad_trials(save_report, chantypes(:,chtype), data_seg);
          all_bad_tri = [all_bad_tri;bt];
      end
      all_bad_tri = unique(all_bad_tri);
      cfg = [];
      cfg.trials = 1:length(data_seg.trial);
      cfg.trials(all_bad_tri) = [];
      cfg.channel = 'meg';
      data_seg = ft_selectdata(cfg,data_seg);
      %% save
      disp(['This data is saved in ' [save_prep subj] ' under the name' subj]);
      ft_write_data([save_prep subj], data_seg, 'dataformat', 'matlab');
  else
      fprintf(fid1,'%s \n',['No MEG data found for subject ' subj]);
      continue
  end
end
