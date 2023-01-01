function [EEG] = rejectbadepochs(EEG, eeg_channel_pos, max_amplitude, min_amplitude)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

epoch_to_reject_ind = zeros(EEG.trials, 1);
for i=1:EEG.trials
    data_i = EEG.data(eeg_channel_pos, :, i);
    if (max(data_i(:)) > max_amplitude) || (min(data_i(:)) < min_amplitude)
        epoch_to_reject_ind(i) = 1;
    end
end
disp(['Found ',num2str(sum(epoch_to_reject_ind)),' epochs to reject']);
if sum(epoch_to_reject_ind) > 0
%     EEG = pop_rejepoch(EEG, epoch_to_reject_ind);
    EEG = pop_rejepoch(EEG, epoch_to_reject_ind, 0);       % added 31-08-2021 to automatically reject epochs without asking for confirmation dialog
end

end

