import os
import mne
import numpy as np
import matplotlib.pyplot as plt
def camcan_preproc(subject):
    data_folder = '/Users/work/camcan'
    data_file = os.path.join(data_folder, subject,
        '_ses-rest_task-rest_proc-raw.fif')
    #data_file = os.path.join(data_folder,
    #                                    'sub-CC _ses-rest_task-rest_proc-raw_sss.fif')
    report = mne.Report(verbose=True)
    report.save('report_basic.html', overwrite=True)
    # read file
    raw = mne.io.read_raw_fif(data_file,preload=True)
    print(raw.info)
    fig = raw.plot_psd(fmax=50,show=False);
    report.add_figs_to_section(fig, captions='Raw',section='PSD')
    meg_channels    = mne.pick_types(raw.info, meg=True, stim=False, ref_meg=False)
    eog1 = mne.pick_channels(raw.ch_names, ['EOG061'])
    eog2 = mne.pick_channels(raw.ch_names, ['EOG062'])
    ecg = mne.pick_channels(raw.ch_names, ['ECG063'])
    mag_channels = mne.pick_types(raw.info, meg='mag')
    # Bandpass between 1 and 45 Hz, bandstop between 50 and 100
    raw.filter(1,45)
    freqs = (50,100)
    raw.notch_filter(freqs=freqs)
    fig = raw.plot_psd(fmax=50,show=False);
    report.add_figs_to_section(fig, captions='After bandpass and notch filters',section='PSD')
    # 2s epochs
    events = mne.make_fixed_length_events(raw,start=0,stop=None,duration=2.0,overlap=0)
    reject = dict(grad=4000e-13, mag=4e-12,eog=250e-6) #from mne tutorial
    epochs = mne.Epochs(raw, events,baseline=None,preload=True, reject=reject)
    # Topoplot before ICA
    fig = epochs.plot_psd_topomap(ch_type='mag', normalize=True,show=False);
    report.add_figs_to_section(fig, captions='Topoplot before ICA',section='Topoplots')
    # ICA to remove artefacts
    from mne.preprocessing import ICA, create_ecg_epochs, create_eog_epochs
    ica = ICA(n_components=0.95, method='fastica').fit(epochs)
    fig = ica.plot_components(show=False);
    ica.exclude = []
    report.add_figs_to_section(fig, captions=('ICA components',' '),section='ICA')
    # Find ECG artefacts
    ecg_epochs = create_ecg_epochs(raw)
    ecg_inds, scores = ica.find_bads_ecg(ecg_epochs, threshold='auto')
    fig = ica.plot_components(ecg_inds,show=False);
    ica.exclude += ecg_inds;
    report.add_figs_to_section(fig, captions='ECG components',section='ICA')
    # Find EOG artefacts
    eog_epochs = create_eog_epochs(raw, tmin=-.5, tmax=.5)
    eog_inds, scores = ica.find_bads_eog(eog_epochs)
    fig = ica.plot_components(eog_inds,show=False);
    ica.exclude += eog_inds;
    report.add_figs_to_section(fig, captions='EOG components',section='ICA')
    # Apply ICA
    cleaned = epochs.copy()
    ica.apply(cleaned)
    fig = cleaned.plot_psd_topomap(ch_type='mag',normalize=True,show=False);
    report.add_figs_to_section(fig, captions='Topoplot after artefact rejection',section='Preprocessed')
    report.save('report.html', overwrite=True)
    # save fif file
    clean_file = os.path.join(data_folder, subject,
        '_cleaned.fif')
    cleaned.save(clean_file, overwrite=True)
