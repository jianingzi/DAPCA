% Calculate the similarity of the test data to the clusters
function [sim] = teWC(Lk, cluster_labels,data, weights,densitythreshold, level, wavename, useSWT)

a=zscore(1:size(data,1))';
num_cells=2*size(data,2);
data=[data,a];
[~,~,~,Lte] = WaveCluster_Preprocess(data, weights, num_cells, densitythreshold, level, wavename, useSWT);

ul= unique(cluster_labels+1);
for i=1:size(ul)
    k=ul(i);
    LK=repmat(Lk(:,:,k,k) , [1,1,num_cells, num_cells]);
    Lte=reshape(Lte,num_cells,[]);LK=reshape(LK,num_cells,[]);
    dot_product = dot(LK, Lte);
    norm_Lk = norm(LK);
    norm_Lte = norm(Lte);
    SIM(i) = norm(dot_product / (norm_Lk * norm_Lte));
end

[~,idx]=max(SIM);
sim=ul(idx)-1;
end