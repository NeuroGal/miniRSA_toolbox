function rdm = calculate_rdm(segs, distance_metric)
% Calculate representational dissimilarity matrix (RDM) of the segments in
% segs (per time-point) based on Matlab's pdist function.
% 
% Example usage: rdm = calculate_rdm(seg, 'correlation')
% Calculates RDM per time-point based on the correlation metric (1 -
% pearson correlation).
%
% Input: segs (dimensions: n_stim * n_features * seg_duration)
%        distance_metric - method to calculate the dissimilarity, e.g.
%        'euclidean', 'correlation', 'cosine', 'mahalanobis' (see options
%        for Matlab's pdist.
% Output: RDMs (dimensions: n_stim * n_trials * seg_duration)
%
% # TODO: 1) add an option for a 'decoding' distance metric (1-accuracy) with
% an option to use decision-value-weighting 2) add an option for cross validated
% distances(see Guggenmos, Sterzer, and Cichy. "Multivariate pattern analysis
% for MEG: A comparison of dissimilarity measures." NeuroImage (2018))
%
% Written by Gal Vishne, lab of Leon Y. Deouell, 2019-2020
% Send bug reports and requests to gal.vishne@gmail.com

cross_validated = false;
if length(size(segs))==4
    cross_validated = true;
    warning('Using cross-validated distances according to dim 4')
end

n_stim = size(segs, 1); seg_duration = size(segs, 3);
rdm = nan(n_stim, n_stim, seg_duration);

if ~cross_validated    
    for t=1:seg_duration
        X = segs(:,:,t);
        rdm(:,:,t) = squareform(pdist(X, distance_metric));
    end
else
    for t=1:seg_duration
        X1 = segs(:,:,t,1);
        X2 = segs(:,:,t,2);
        rdm(:,:,t) = (pdist2(X1,X2, distance_metric)+pdist2(X2,X1, distance_metric))/2;
    end   
end
    
end