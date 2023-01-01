function [] = erp_subplot(EEG_epoch, chan_num, plot_ci, col)
%[] = erp_subplot(EEG_epoch, chan_num, plot_ci, col)
%   Function used to plot ERP

t = EEG_epoch.times(:);
data_epoch  = squeeze(EEG_epoch.data(chan_num,:,:));
data_mean   = mean(data_epoch, 2);
data_std    = std(data_epoch, 0, 2);
n_trials    = EEG_epoch.trials;
CI_upper    = data_mean + ((1.96.*data_std./sqrt(n_trials)));
CI_lower    = data_mean - ((1.96.*data_std./sqrt(n_trials)));

% Plot ERP 
plot(t, data_mean, 'Color', col, 'LineWidth', 2); hold on;
% Plot Confidence Interval
if plot_ci
    fill([t; flipud(t)], [CI_upper; flipud(CI_lower)], col, 'EdgeColor', 'none','FaceAlpha', 0.2);
end
xlabel('Time (ms)');
ylabel('Amplitude (uV)');

end

