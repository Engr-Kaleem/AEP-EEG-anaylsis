function [ freq_band_power ] = freqanalysis_computebandenergy(pxx, f_vect, freq_bands)
%[ freq_band_power ] = freqanalysis_computebandenergy(pxx, f_vect, freq_bands)
%   Compute the mean power by frequency bands

n_bands = length(freq_bands);
freq_band_power = zeros(n_bands, 1);
for i = 1:n_bands
    f_band_i = freq_bands{i};
    if size(f_band_i, 1) == 1
        f_band_vect = (f_vect > f_band_i(1)) & (f_vect < f_band_i(2));
    else
        f_band_vect =  (f_vect > f_band_i(1, 1)) & (f_vect < f_band_i(1, 2));
        for j = 2:size(f_band_i, 1)
            f_band_vect = f_band_vect | ((f_vect > f_band_i(j, 1)) & (f_vect < f_band_i(j, 2)));
        end
    end
    freq_band_power(i) = mean(pxx(f_band_vect));   
end

end

