% % The imported dataset "dataset.mat" is obtained by adding a -10% drift fault to the dataset.

% % The performance between different methods needs to be compared under the same drift fault dataset.

% % The following is the set of parameters to be debugged in the experiment.
% % alpha------------A percentile for calculating the density threshold
% % omege------------The smoothing window size
% % beta------------- The confidence level for determining the detection threshold

clc
clear
close all
r_values = -10;
    for i= 1:length(r_values)      
        disp(['------------------------- drift fault rate is ' num2str(r_values(i)) '%. -------------------------'])
        % load dataset
        load dataset.mat    % This dataset contains the following data:
        % 1. A vector of label columns named  "label" (containing 0 and 1, with 0 being normal and 1 being faulty);
        % 2. m columns of the matrix named "normal_data1" represent normal training data where m represents the number of sensors;
        % 3. m columns of the matrix named "uncertain_data1" represents the testing data to be detected where m represents the number of sensors.

        % parameter setting  
        alpha=0.95*100; % for calculating the density threshold
        omega=15;   % the smoothing window size 
        beta=0.99;  % he confidence level

        % data pre-processing
        mu_train = mean(normal_data1, 1);
        sigma_train = std(normal_data1, 0, 1);
        xtrain=(normal_data1 - mu_train) ./ sigma_train;
        xtest = (uncertain_data1 - mu_train) ./ sigma_train;

        % temporal waveclusting algorithm for data clusting
        [cluster_labels,~,~,~,~,sigcells] = WaveCluster(xtrain, [], 'alpha%', 1, 'bior2.2', 1);
        similarity=teWC(sigcells,cluster_labels,xtest, [],'alpha%', 1, 'bior2.2', 1);
        indices = find(similarity == cluster_labels);
        Xtrain = xtrain(indices, :);

        % a roubust PCA-based drift detection method
        [FAR(i),F1(i),fai{i},kesi{i},result{i},AUC(i)]=RPCA(Xtrain,xtest,label,omega,beta,i);

        % output the result
        disp(['FAR = ' num2str(FAR(i), '%.2f') '% , F1 = ' num2str(F1(i), '%.2f') '% , and AUC = ' num2str(AUC(i),'%.4f')]);
        subplot(length(r_values),1,i)
        plot((1:size(xtest)),uncertain_data1(:,2),'b-*')
        hold on
        ikk=find(result{i}==1);
        scatter(ikk,uncertain_data1(ikk,2),30,'ro','filled');
        hold on
        plot([find(label==1,1),find(label==1,1)],[min(uncertain_data1(:,2)),max(uncertain_data1(:,2))],'k--');
        hold on
        legend('Sensor output','Detected faulty','Fault start point','Location','best','NumColumns',3);
        xlabel('Testing sample');
        ylabel('Sensor output');
        title(['Drift rate: ' num2str(r_values(i)) '%']);
    end
