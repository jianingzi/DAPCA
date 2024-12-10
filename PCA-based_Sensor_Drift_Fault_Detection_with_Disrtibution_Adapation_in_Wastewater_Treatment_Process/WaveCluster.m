% Run the WaveCluster algorithm on the input data. Returns the cluster labels.
% Accepts the set of observations in the rows and features in the columns (including x,y coordinates), an optional vector of observation weights, the number of cells to be divided, the
% count threshold above which the cell is considered ‘significant’, the level of wavelet decomposition to be analysed, and the name of the wavelet to be used. The name of the wavelet to be used.
% If useSWT is non-zero, the stationary wavelet transform will be used instead of the discrete wavelet transform. This improves performance for small data sets.
% The threshold, level and wavelet name are optional; the default values are 10% of the data set, 1 and 2,2 for bi-orthogonal wavelets.
% The cell size can be a scalar or a vector of the same size as the number of features in the data matrix. The number of features in the data matrix is the same.

function [cluster_labels, sigcells,clustergrid, counts, datacellindices, wdata] = WaveCluster(data, weights, densitythreshold, level, wavename, useSWT)
a=zscore(1:size(data,1))';
num_cells=2*size(data,2);
data=[data,a];
[sigcells, datacellindices, counts, wdata] = WaveCluster_Preprocess(data, weights, num_cells, densitythreshold, level, wavename, useSWT);
clustergrid = bwlabeln(sigcells);
linidx = num2cell(datacellindices, 1);
cluster_labels=clustergrid(sub2ind(size(clustergrid), linidx{:}));
%     t=0;
%     t_label=[];
%     for i=1:max;
%         ismember(cluster_labels,i)
%         t=t+1;
%         t_label=[t_label,i];
%     end
%
%  for j=1:siez(t_label,2)
%  for i=1:size(cluster_labels,1)
%      t_label(:,j)=
%  disp(['Discovered ' num2str(max(cluster_labels)) ' clusters.']);
end