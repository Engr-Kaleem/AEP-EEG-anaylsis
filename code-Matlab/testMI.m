% Load Fisher's iris dataset (meas = samples, species = class labels)
close all;
clear all;
load fisheriris meas species 

% Convert class label strings into integer labels
specs = unique(species);
labels = zeros(size(species));
for k = 1:length(specs)
    labels(ismember(species,specs(k))) = k;
end

% Add Gaussian noise to the measurement data (original fisheriris is too easy for classification).
noiselevel = 1;
meas = meas+randn(size(meas)).*noiselevel;

% Generate a set of new features through random projection from the
% original 4 features to d dimensions.
d = 200;
M = randn(size(meas,2),d);
M = sqrt(ones./(sum((M.*M)')))'*ones(1,size(M,2)).*M; % Normalize M rows
features = meas*M;

% Replace max 50% of the generated features with random noise features
a = 1 + floor(size(features,2)*rand(round(d/2),1));
features(:,a) = randn(size(features,1),length(a));


fprintf('Feature selection using SD, MI, RSFS, SFS and SFFS\n');
fprintf('Please see the source code for more information\n');
fprintf('Evaluation started\n');

echo on

%% Select features using different algorithms
[F_MI,W_MI] = MI(features,labels,3);
