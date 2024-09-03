function resample_data = Func_Resample(original_data, origianl_fs, CutOff_fz)
    % Original signal and sampling rate
    original_signal = original_data;
    original_sampling_rate = origianl_fs; % Hz
    % Define the desired frequency cutoff
    desired_cutoff_frequency = CutOff_fz; % Hz
    % Determine the downsampling factor
    downsampling_factor = original_sampling_rate / (2 * desired_cutoff_frequency);
    % Downsample the signal
    resample_data = downsample(original_signal, downsampling_factor);
    % Determine the new sampling rate
    new_sampling_rate = original_sampling_rate / downsampling_factor;
    % Display the new sampling rate
    disp(['New sampling rate: ', num2str(new_sampling_rate), ' Hz']);
end