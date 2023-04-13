function [features,weights] = MI(features,labels,Q)
% function [features,weights] = MI(features,labels,Q)
% Estimates the mutual information between features and associated class labels using a quantized feature space.
%
%   Inputs:
%           features:       N x F sized matrix of features, where N is the number of samples and F is the number of features        
%           labels:         N x 1 sized vector of class labels corresponding to each sample
%           Q:              the number of quantization levels used for the features (default = 12)
%
%   Outputs:
%           features:       F x 1 sized vector of feature indices in the
%                              descending order of relevance.
%           weights:        F x 1 sized vector of feature relevances (MIs) in the
%                              descending order.
%
% Author: Okko Rasanen, 2013. Mail: okko.rasanen@aalto.fi
%
% The algorithm can be freely used for research purposes. 
%
% Please see J. Pohjalainen, O. Rasanen & S. Kadioglu: "Feature Selection Methods and 
% Their Combinations in High-Dimensional Classification of Speaker Likability, 
% Intelligibility and Personality Traits", Computer Speech and Language, 2015, for more details.

if nargin <3
    Q = 12;
end

edges = zeros(size(features,2),Q+1);

% Compute feature-specific quantization bins so that each bin has approximately equal number of
% samples in the training set
for k = 1:size(features,2)
    
    minval = min(features(:,k));
    maxval = max(features(:,k));
    if minval==maxval
        continue;
    end
    
    quantlevels = minval:(maxval-minval)/500:maxval;
    
    N = histc(features(:,k),quantlevels);
    
    totsamples = size(features,1);
    
    N_cum = cumsum(N);
    
    edges(k,1) = -Inf;
    
    stepsize = totsamples/Q;
    
    for j = 1:Q-1
        a = find(N_cum > j.*stepsize,1);
        edges(k,j+1) = quantlevels(a);
    end
    
    edges(k,j+2) = Inf;
end

% Quantize data according to the obtained bins
S = zeros(size(features));
for k = 1:size(S,2)
    S(:,k) = quantizedvec(features(:,k),edges(k,:))+1;   
end

% Compute mutual information (MI) between the quantized features and
% the class labels
I = zeros(size(features,2),1);
for k = 1:size(features,2)   
    I(k) = computeMI(S(:,k),labels,0);
end

% Sort features into descending order

[weights,features] = sort(I,'descend');