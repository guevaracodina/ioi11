function yfit = ...
    crossfun(xtrain,ytrain,xtest,rbf_sigma,boxconstraint)
% Write a function called crossfun to calculate the predicted
% classification yfit from a test vector xtest, when the SVM is trained on
% a sample xtrain that has classification ytrain. Since you want to find
% the best parameters rbf_sigma and boxconstraint, include those in the
% function.

% Train the model on xtrain, ytrain, 
% and get predictions of class of xtest
svmStruct = svmtrain(xtrain,ytrain,'Kernel_Function','rbf',...
   'rbf_sigma',rbf_sigma,'boxconstraint',boxconstraint);
yfit = svmclassify(svmStruct,xtest);

end

