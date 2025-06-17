function [kfoldLoss, predictions] = lda_kfold_cv(X, Y, gamma,k,desiredCPtoKeep)
% Performs k-fold cross-validation using Linear Discriminant Analysis (LDA)
% 
% Inputs:
%   X - Feature matrix (samples x features)
%   Y - Target labels (samples x 1)
%   k - Number of folds
%
% Outputs:
%   accuracy - Average classification accuracy across folds
%   predictions - Predicted labels for each fold (in original order)

    % Ensure Y is a column vector
    Y = Y(:);
    N = size(X, 1);

    % Create cross-validation partition
    cv = cvpartition(Y, 'KFold', k);

    predictions = nan(N, 1);
    acc = zeros(k, 1);

    for i = 1:k
        trainIdx = training(cv, i);
        testIdx = test(cv, i);

        % Train LDA
        X0=X(:,trainIdx);
        badRows=max(X0,[],2)==0;
        if sum(badRows)<=(size(X0,1)-desiredCPtoKeep)
        model = fitcdiscr(X0(~badRows,:)', Y(trainIdx),'DiscrimType','linear','Gamma',gamma);
        Ypred = predict(model, X(~badRows,testIdx)');
        else
        model = fitcdiscr(X0', Y(trainIdx),'DiscrimType','diaglinear');
        Ypred = predict(model, X(:,testIdx)');
        end
        % Predict
        % Ypred = predict(model, X(~badRows,testIdx)');
        predictions(testIdx) = Ypred;

        % Accuracy
        acc(i) = mean(Ypred == Y(testIdx));
    end

    % Average accuracy
    kfoldLoss = 1-mean(acc);
end