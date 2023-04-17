close all;
clear all;

load features.mat

stftfeat=zeros(sum(epochs),size(smat,1)*size(smat,2));
erspfeat=zeros(sum(tlabels),size(erspmat,1)*size(erspmat,2));
itcfeat=zeros(sum(tlabels),size(itcmat,1)*size(itcmat,2));

ind=1;
for i=1:length(epochs);
   for ep=1:epochs(i);
      s=smat(:,:,ep,i) ;   
    
      stftfeat(ind,:)= reshape(s',1,numel(s));
      stft_label(ind)=labels_stft(i,ep);
      ind=ind+1;
   end

end


ind=1;
for i=1:length(tlabels);
   for ep=1:tlabels(i);
      
      e=erspmat(:,:,ep,i); 
      it=itcmat(:,:,ep,i);
      erspfeat(ind,:)= reshape(e',1,numel(e));
      itcfeat(ind,:)= reshape(it',1,numel(it));
      all_label(ind)=labels(i,ep);
      ind=ind+1;
   end

end



%% Compute mutual information between each feature set and the target

stft_n=normalize(stftfeat);
ersp_n=normalize(erspfeat);
itc_n=normalize(itcfeat);



labels=all_label';
for k = 1:size(ersp_n,2)   
    erspMI(k) = computeMI(ersp_n(:,k),all_label');
end

for k = 1:size(itc_n,2)   
    itcMI(k) = computeMI(itc_n(:,k),all_label');
end

for k = 1:size(stft_n,2)   
    stftMI(k) = computeMI(stft_n(:,k),stft_label');
end
%%

disp(['MI ersp: ', num2str(abs(sum(erspMI))/length(all_label))]);
disp(['MI ITC: ', num2str(abs(sum(itcMI))/length(all_label))]);
disp(['MI STFT: ', num2str(abs(sum(stftMI))/length(stft_label))]);
%%
% 
% 
% 
% MI(features,labels,3)
% 
% %%
% mi = zeros(1, 3);
% for i = 1:3
%     % Compute joint histogram of feature set i and target
%     edgesX = linspace(min(X{i}(:)), max(X{i}(:)), 10); % edges of histogram bins for feature set i
%     edgesY = unique(Y); % unique values of target
%     hyx = histcounts2(Y, X{i}, edgesY, edgesX); % joint histogram of Y and feature set i
%     hy = histcounts(Y, edgesY); % histogram of Y
%     hx = histcounts(X{i}, edgesX); % histogram of feature set i
%     pxy = hyx / sum(hyx(:)); % joint probability distribution of Y and feature set i
%     py = hy / sum(hy); % marginal probability distribution of Y
%     px = hx / sum(hx); % marginal probability distribution of feature set i
%     mi(i) = sum(sum(pxy .* log2(pxy ./ (px' * py)))); % mutual information between Y and feature set i
% end
% 
% % Select feature set with highest mutual information
% [~, idx] = max(mi); % index of feature set with highest mutual information
% if idx == 1
%     selected_features = X1;
% elseif idx == 2
%     selected_features = X2;
% elseif idx == 3
%     selected_features = X3;
% end
% 
% 

%%
save('featuremat.mat', 'stftfeat', 'erspfeat','itcfeat','all_label','epochs','stft_label','tlabels'); % save both variables to a file