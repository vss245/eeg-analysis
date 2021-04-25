% EEG classification of motor imagery task
%% Loading the data
load('mi_data.mat');
eeg_data=fv_imag.x;
ntrials = size(eeg_data,3)/2;
nchan = size(eeg_data,2);
nsam = size(eeg_data,1);
fs=100;
%% Filtering
% known changes in activity occur in the mu band (8-13 Hz)
% band-pass filter the data in this range
nfs = fs/2;
[b, a] = butter(4, [8/nfs 13/nfs]);
feeg_data = filtfilt(b,a,eeg_data);
% split the filtered data into training and testing sets
train_set= feeg_data(:,:,1:75);
test_set=feeg_data(:,:,76:150);
%% Step 2 - Calculating the CSP
% get hand labels and split the training data into left and right hand
h_labels = (fv_imag.y(1,1:75))'; % 1 for left hand, 0 for right hand
lh_data = train_set(:,:,find(h_labels));
rh_data = train_set(:,:,find(~h_labels));
% obtain channel covariance matrices (A'A) for each trial and average 
% prealloc
r_tri = size(rh_data,3);
l_tri = 75-r_tri;
rh_cov = ones(nchan,nchan,r_tri);
lh_cov = ones(nchan,nchan,l_tri);
% right hand
for i = 1:r_tri
    trial = rh_data(:,:,i);
    % subtract the mean
    trial = bsxfun(@minus,trial,squeeze(mean(trial,1)));
    cov_mat = cov(trial)/trace(cov(trial)); 
    rh_cov(:,:,i) = cov_mat;
end
C1 = mean(rh_cov,3); % average over trials
% left hand
for i = 1:l_tri
    trial = lh_data(:,:,i);
    trial = bsxfun(@minus,trial,squeeze(mean(trial,1)));
    cov_mat = cov(trial)/trace(cov(trial));
    lh_cov(:,:,i) = cov_mat;
end
C2 = mean(lh_cov,3); % average over trials 
% generalized eigenvalue problem: C1*w = lambda*C2*w
[ei_vec, ei_val] = eig(C1, C1+C2);
% sort eigenvalues
[ei_val_d,indexes] = sort(diag(ei_val),'descend');
ei_vec = ei_vec(:, indexes);
%retaining 6 spatial filters (3 with highest eigenvalues, 3 with lowest)
csp = [ei_vec(:,1:3) ei_vec(:,82:84)];
%% Step 3 - Project the training data with the CSP filters
proj_mat = ones(nsam,6,ntrials);
for n = 1:ntrials
   proj_mat(:,:,n) = train_set(:,:,n)*csp;
end
%variance across time and log
logvar = squeeze(log(var(proj_mat)));
%% Step 4 - LDA
%train LDA
%split by class and calculate means
lv_left = logvar(:,find(h_labels));
lv_right = logvar(:,find(~h_labels));
mu1 = mean(lv_left,2);
mu2 = mean(lv_right,2);
w = cov(logvar')\(mu1-mu2); %LDA weights
b = (mu1+mu2)'*w/2; %bias
%% Step 5 - Apply LDA to training data
tr_lda = logvar'*w-b;
scores = sign(tr_lda);
true_lab = h_labels;
true_lab(true_lab==0)=-1;
%accuracy for the training set
acc = (sum(scores==true_lab)/length(scores));
disp(['Training accuracy: ', num2str(acc*100),'%'])
%% Step 6 - Apply CSP and LDA to testing data
proj_mat_test = ones(nsam,6,ntrials);
for n = 1:ntrials
   proj_mat_test(:,:,n) = test_set(:,:,n)*csp;
end
logvar_test = squeeze(log(var(proj_mat_test)));
test_lda = logvar_test'*w-b;
scores_test = sign(test_lda);
true_lab_test = fv_imag.y(1,76:150)';
true_lab_test(true_lab_test==0)=-1;
%accuracy for the testing set
acc_test = (sum(scores_test==true_lab_test)/length(scores));
disp(['Testing accuracy: ', num2str(acc_test*100),'%'])
%% Step 7 - Plotting based on scalp topographies
%chanlocs = struct('labels',fv_imag.clab);
%chanlocs = pop_chanedit(chanlocs); %get standard channel locations from
%eeglab using MNI coordinates
load('chanlocs'); % 84-electrode channel locations from eeglab
for n = 1:ntrials
    trial = feeg_data(:,:,n);
    trial = bsxfun(@minus,trial,squeeze(mean(trial,1)));
    cov_mat = cov(trial)/trace(cov(trial));
    all_cov(:,:,i) = cov_mat;
end
sigma_x = mean(all_cov,3);
sigma_s = cov(logvar');
patterns = (sigma_x*csp)/sigma_s;
figure
for i = 1:6
    subplot(2,3,i)
    topoplot(patterns(:,i),chanlocs,'electrodes','on');
    title(i);
end
sgtitle('Activation patterns')
