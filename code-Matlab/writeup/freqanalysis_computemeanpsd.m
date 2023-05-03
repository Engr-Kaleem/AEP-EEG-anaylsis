function [FFT_i,pxx_db_mean, f_vect] = freqanalysis_computemeanpsd(EEG, chan_num, nfft, t_min, t_max)
%[pxx_db_mean, f_vect] = freqanalysis_computemeanpsd(EEG, chan_num, nfft, t_min, t_max)
%  Compute the power spectral density for each trial and compute the mean 
%  The PSD is estimated using fft 
%  Use a hanning window

data = squeeze(EEG.data(chan_num, :, :));
t_ind       = (EEG.times > 1000*t_min) & (EEG.times < 1000*t_max);
window_fft  = hanning(sum(t_ind)); 
% window_fft  = ones(sum(t_ind),1);
nfft        = 4096;
pxx_db_all  = zeros(1+nfft/2, EEG.trials);
% For each trial
for i = 1:EEG.trials
    trial_i     = data(t_ind, i);
    fft_i       = fft(trial_i.*window_fft, nfft);
    pxx_db_all(:, i)  = 10*log10(abs(fft_i(1:1+nfft/2)).^2);
    FFT_i(:,i)=fft_i;
end
pxx_db_mean = mean(pxx_db_all, 2);
f_vect = linspace(0, EEG.srate/2, 1+nfft/2);
end

