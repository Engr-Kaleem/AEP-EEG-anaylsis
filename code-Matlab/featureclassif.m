close all;
clear all;
load featuremat.mat

stft_n=(stftfeat);
ersp_n=(erspfeat);
itc_n=(itcfeat);
Y=all_label';

%% ersp  feature for classification
X=ersp_n;


% Split the dataset into training and testing sets
cv = cvpartition(size(X,1),'HoldOut',0.2);
X_train = X(training(cv),:);
Y_train = Y(training(cv));
X_test = X(test(cv),:);
Y_test = Y(test(cv));

% SVM classification
svm_model = fitcecoc(X_train, Y_train);
svm_predicted_labels = predict(svm_model, X_test);



%% LDA classification
lda_model = fitcdiscr(X_train, Y_train);
lda_predicted_labels = predict(lda_model, X_test);



%% Compute the classification accuracy of each model
svm_accuracy = sum(svm_predicted_labels == Y_test)/length(Y_test);
lda_accuracy = sum(lda_predicted_labels == Y_test)/length(Y_test);

% % SVM confusion matrix
svm_cmat_ersp = confusionmat(Y_test, svm_predicted_labels);
% figure();
% heatmap(svm_cmat, {'inphase ', 'antiphase'}, {'inphase', 'antiphase'});
% title('Confusion Matrix for SVM Classifier ERSP');
% xlabel('Predicted Labels');
% ylabel('True Labels');
% 
% % LDA confusion matrix
lda_cmat_ersp = confusionmat(Y_test, lda_predicted_labels);
% figure();
% heatmap(lda_cmat, {'inphase ', 'antiphase'}, {'inphase', 'antiphase'});
% title('Confusion Matrix for LDA Classifier ERSP');
% xlabel('Predicted Labels');
% ylabel('True Labels');


disp(['SVM accuracy ersp: ', num2str(svm_accuracy)]);
% disp(['k-NN accuracy: ', num2str(knn_accuracy)]);
disp(['LDA accuracy ersp: ', num2str(lda_accuracy)]);

%%
X=itc_n;


% Split the dataset into training and testing sets
cv = cvpartition(size(X,1),'HoldOut',0.2);
X_train = X(training(cv),:);
Y_train = Y(training(cv));
X_test = X(test(cv),:);
Y_test = Y(test(cv));

% SVM classification
svm_model = fitcecoc(X_train, Y_train);
svm_predicted_labels = predict(svm_model, X_test);



% LDA classification
lda_model = fitcdiscr(X_train, Y_train);
lda_predicted_labels = predict(lda_model, X_test);


%%
% Compute the classification accuracy of each model
svm_accuracy = sum(svm_predicted_labels == Y_test)/length(Y_test);
lda_accuracy = sum(lda_predicted_labels == Y_test)/length(Y_test);



% % SVM confusion matrix
svm_cmat_itc = confusionmat(Y_test, svm_predicted_labels);
% figure();
% heatmap(svm_cmat, {'inphase ', 'antiphase'}, {'inphase', 'antiphase'});
% title('Confusion Matrix for SVM Classifier ITC');
% xlabel('Predicted Labels');
% ylabel('True Labels');
% 
% % LDA confusion matrix
lda_cmat_itc = confusionmat(Y_test, lda_predicted_labels);
% figure();
% heatmap(lda_cmat, {'inphase ', 'antiphase'}, {'inphase', 'antiphase'});
% title('Confusion Matrix for LDA Classifier ITC');
% xlabel('Predicted Labels');
% ylabel('True Labels');

disp(['SVM accuracy itc: ', num2str(svm_accuracy)]);
disp(['LDA accuracy itc: ', num2str(lda_accuracy)]);
%%

%%
X=stft_n;
Y=stft_label';

% Split the dataset into training and testing sets
cv = cvpartition(size(X,1),'HoldOut',0.2);
X_train = X(training(cv),:);
Y_train = Y(training(cv));
X_test = X(test(cv),:);
Y_test = Y(test(cv));




% SVM classification
svm_model = fitcsvm(X_train,Y_train,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto');;
svm_predicted_labels = predict(svm_model, X_test);



% LDA classification
lda_model = fitcdiscr(X_train, Y_train);
lda_predicted_labels = predict(lda_model, X_test);



% Compute the classification accuracy of each model
svm_accuracy = sum(svm_predicted_labels == Y_test)/length(Y_test);
lda_accuracy = sum(lda_predicted_labels == Y_test)/length(Y_test);

% % SVM confusion matrix
svm_cmat_stft = confusionmat(Y_test, svm_predicted_labels);
% figure();
% heatmap(svm_cmat, {'inphase ', 'antiphase'}, {'inphase', 'antiphase'});
% title('Confusion Matrix for SVM Classifier  STFT');
% xlabel('Predicted Labels');
% ylabel('True Labels');
% 
% % LDA confusion matrix
lda_cmat_stft = confusionmat(Y_test, lda_predicted_labels);
% figure();
% heatmap(lda_cmat, {'inphase ', 'antiphase'}, {'inphase', 'antiphase'});
% title('Confusion Matrix for LDA Classifier STFT');
% xlabel('Predicted Labels');
% ylabel('True Labels');

disp(['SVM accuracy stft: ', num2str(svm_accuracy)]);

disp(['LDA accuracy stft: ', num2str(lda_accuracy)]);
%%
% data=ersp_n;
% labels=all_label';
% 
% 
% % Split the data into training and testing sets
% train_size = 0.7; % percentage of data used for training
% rng(1); % set the random seed for reproducibility
% idx = randperm(size(data,1));
% train_data = data(idx(1:round(train_size*size(data,1))),:);
% train_labels = labels(idx(1:round(train_size*size(data,1))),:);
% test_data = data(idx(round(train_size*size(data,1))+1:end),:);
% test_labels = labels(idx(round(train_size*size(data,1))+1:end),:);
% 
% 
% num_trees = 50;
% model = TreeBagger(num_trees, train_data, train_labels);
% 
% % Predict on the testing set
% predictions = predict(model, test_data);
% 
% % Evaluate the performance
% accuracy = sum(str2double(predictions)==test_labels)/length(test_labels);
% disp(['Accuracy: ' num2str(accuracy)]);