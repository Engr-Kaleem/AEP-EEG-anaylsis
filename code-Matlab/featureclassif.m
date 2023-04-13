close all;
clear all;
load featuremat.mat

stft_n=normalize(stftfeat);
ersp_n=normalize(erspfeat);
itc_n=normalize(itcfeat);
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



% Compute the classification accuracy of each model
svm_accuracy = sum(svm_predicted_labels == Y_test)/length(Y_test);
lda_accuracy = sum(lda_predicted_labels == Y_test)/length(Y_test);

disp(['SVM accuracy itc: ', num2str(svm_accuracy)]);
disp(['LDA accuracy itc: ', num2str(lda_accuracy)]);
%%

%%
X=stft_n;


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



% Compute the classification accuracy of each model
svm_accuracy = sum(svm_predicted_labels == Y_test)/length(Y_test);
lda_accuracy = sum(lda_predicted_labels == Y_test)/length(Y_test);

disp(['SVM accuracy stft: ', num2str(svm_accuracy)]);

disp(['LDA accuracy stft: ', num2str(lda_accuracy)]);
%


