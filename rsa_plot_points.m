function fig = rsa_plot_points(segs, method, time_ms, varargin)
% Plot a representation of the given segments in 2D at specific
% time-points given by time_ms. Each time-point is plotted in a different
% subplot, where each segment is represented as a point.
%
% Default settings for conversion of ms to indices in data:
%   sampling_rate   = 250 Hz
%   baseline        = 100 ms
% Default distance metric for 'tsne', 'mds' and 'cmds' is 'euclidean'.
% (to change see optional input)
%
% Example usage (more options under input section):
%   - fig = rsa_plot_points(segs, 'tsne', [-50:50:300], 'link_plots')
%       Plot segments using tsne algorithm with euclidean distance (default)
%       at time-points -50ms to 300ms. Conversion to indices is done using
%       default settings. Feed each time-point as an initial estimate for
%       the algorithm for the next point and synchronize axes limits.
%   - fig = rsa_plot_points(segs, 'pca' , [-50:50:300], 'categories', cat_inds, cat_names)
%       Use pca to plot the points, and color-code according to the
%       categories in cat_inds.
%   - fig = rsa_plot_points(segs, 'mds', [], 'baseline', 200, 'metric', 'cosine')
%       Plot segments using non-classical mds with cosine dissimilarity at
%       the default time-points (see below), and assume baseline is 200ms
%       for conversion to indices.
% 
% Input:
%   segs (segments to plot, dimensions: n_stim * n_features * seg_duration)
%   method - how to perform dimensionality reduction, options: 'tsne',
%       'pca', 'cmds' (classical mds), 'mds' (non-classical mds)
%       * 'tsne', 'mds', 'cmds' are performed using euclidean distance
%         unless provided with a 'metric' optional argument.
%       * tsne may take time to compute - therefore the code provides
%         feedback regarding the current plot in that case.
%   time_ms (optional): plot the RDMs corresponding to these time-points
%       (in ms around stimulus onset, e.g. 50ms before onset is -50, and
%       50ms after stimulus onset is +50).
%       Default conversion to indices in rdm is done using sampling rate
%       250 Hz and 100ms baseline. Change using optional arguments.
%       If no time_ms is supplied the default is 6 evently spaced time-points
%       where the distance from segment boundaries is 1/2 distance between
%       time-points, e.g. if the segment is 300 samples the default is to
%       plot samples 25, 75, 125, 175, 225 and 275. 
%   Additional optional arguments (varargin):
%       - 'metric' followed by a name of a distance metric. Method to
%           calculate distances for 'tsne', 'mds', 'cmds'. Recommended to use
%           the same metric used for calculate_rdm. Examples: 'euclidean',
%           'correlation', 'cosine', 'mahalanobis' (for more options see 
%           Matlab's pdist). Default is euclidean.
%       - 'n_pc' followed by integer - number of pcs to use for 'tsne'.
%           Default is min(50, n_features)
%       - 'color_code' followed by a numeric array size n_stim * 1\2 (!) -
%           color code for the points denoting the different segments
%           (shared across all segments). column 1 - color of the points,
%           column 2 (if supplied) - color for surrounding circles (e.g.
%           use the 1st color scheme to mark categories and the 2nd results
%           of a clustering algorithm).
%           Both column values are used as indices to jet color map. 
%           Default color code: each stimulus has it's own color (value are
%           1 to n_stim in the jet color map)
%       - 'link_plots' - if method = 'pca', perform pca on all
%           time-points together (default is separately for each point).
%           If method = 'mds' or 'tsne' the previous time-point solution
%           is fed as an initial estimate to the next point. This eases
%           comparison between time-points, but isn't recommended in case
%           of large temporal gaps or fast changes since the algorithm may
%           not converge.
%           In all cases link_plots synchronizes the axes limits for all
%           time-points.
%       - 'categories' - color the segments according to categories and add
%           labels. 'categories' must be followed by (in this order!):
%           1) cat_inds - a numeric vector (length n_stim) indicating
%              the category assignment of each stimulus. Use numbers from 1
%              to n_categories.
%           2) cat_names - a cell array of the name of each category.
%           In case color_code is also given cat_inds are used for the
%           points and color_code for the surrounding circles. 
%       - 'sampling_rate' followed by a number (in Hz) - change the
%           sampling rate for the conversion of ms to indices.
%       - 'baseline' followed by a number (in ms) - change the baseline
%           seg_duration for the conversion of ms to indices.
%
% Output: figure with 2D representation of the given segments at specific
%   time-points. Each time-point is plotted in a different subplot, where
%   each segment is represented as a point.
%
% Written by Gal Vishne, lab of Leon Y. Deouell, 2019-2020
% Send bug reports and requests to gal.vishne@gmail.com

n_stim = size(segs,1); n_features = size(segs,2); seg_duration = size(segs,3);

link_plots = false;
categories = false;
baseline = 100;         % in ms
sampling_rate = 250;    % in Hz
n_pc = min(n_features,50);
color_code = (1:n_stim)';
sub_color_code = false;
distance_metric = 'euclidean';
arg  = 1;
while arg <= size(varargin,2)
    switch varargin{arg}
        case 'categories'
            categories = true;
            cat_inds = varargin{arg+1};
            if length(cat_inds) ~= n_stim
                error('Category indices vector should match the number of stimuli');
            end
            if size(cat_inds,2)>size(cat_inds,1)
                cat_inds = cat_inds';
            end
            cat_names = varargin{arg+2};
            if ~iscell(cat_names)
                error('Category names should be given in a cell array')
            end
            if size(cat_names,2)>size(cat_names,1)
                cat_names = cat_names';
            end
            if length(unique(cat_inds))~=length(cat_names)
                error("Number of category names doesn't match number of categories")
            end
            arg = arg + 3;
        case 'baseline'
            baseline = varargin{arg+1};
            arg = arg + 2;
        case 'sampling_rate'
            sampling_rate = varargin{arg+1};
            arg = arg + 2;
        case 'color_code'
            color_code = varargin{arg+1};
            if size(color_code,2)>2
                error('No more than 2 color-codes')
            elseif size(color_code,1)~=n_stim
                error('Color-code length should match the number of stimuli')
            end
            sub_color_code = true;
            arg = arg + 2;
        case 'link_plots'
            link_plots = true;
            arg = arg + 1;
        case 'metric'
            distance_metric = varargin{arg+1};
            arg = arg + 2;
        case 'n_pc'
            n_pc = varargin{arg+1};
            if n_features < n_pc
                warning('Using %d pc for tsne (n_pc needs to be smaller than the number of features)', n_features)
                n_pc = n_features;
            end
            arg = arg + 2;
        otherwise
            error(['Unknown optional argument name: ' varargin{arg} '.']);
    end
end

if exist('time_ms','var')
    if isempty(time_ms)
        time_inds = linspace(1, seg_duration, 13);
        time_inds = round(time_inds(2:2:13));
        time_ms = ((time_inds-1)*1000/sampling_rate - baseline);
    else
        time_inds = floor((time_ms+baseline)*sampling_rate/1000)+1; %turn ts into indices
        if min(time_inds) < 1 || max(time_inds) > seg_duration
            error('Given temporal range exceeds segment duration');
        end
    end
else
    time_inds = linspace(0, seg_duration, 13);
    time_inds = round(time_inds(2:2:13));
    time_ms = ((time_inds-1)*1000/sampling_rate - baseline);
end

if length(time_inds)<3
    num_rows=1;
else
    num_rows=2;
end
num_cols=round(length(time_inds)/num_rows);
fig = figure('Units','Normalized','Position',[0.2 0.3 0.2+0.07*num_cols 0.4]);

if categories
    if size(color_code,2)==2
        warning('Using only the 1st color code (basic color-code is by categories)')
    end
    if sub_color_code 
        color_code(:,2) = color_code(:,1);
    end
    color_code(:,1) = cat_inds;
end

if strcmp(method, 'tsne')
    fprintf('starting tsne calculations:\n')
elseif strcmp(method, 'pca') && link_plots
    dataforpca = reshape(permute(segs(:,:,time_inds),[2,1,3]),n_features,n_stim*length(time_inds))';
    [~,scores] = pca(dataforpca, 'NumComponents', 2);
end

ha = [];
for t=1:length(time_inds)
    w_inset = [0.1 0.02]; h_inset = [0.07 0.07]; w_gap = 0.018; h_gap = 0.08;
    width = (1-sum(w_inset)-(num_cols-1)*w_gap)/num_cols; height = (1-sum(h_inset)-h_gap)/num_rows;
    if t<=num_cols
        pos = [w_inset(1)+(t-1)*(width + w_gap) h_inset(1)+(h_gap+height)*(num_rows-1) width height];
    else
        pos = [w_inset(1)+(t-1-num_cols)*(width + w_gap) h_inset(1) width height];
    end
    ha(t) = axes('Position',pos);
    
    t0=time_inds(t);
    data_toplot = segs(:,:,t0);
    if strcmp(method, 'mds') || strcmp(method, 'cmds')
        data_toplot = squareform(pdist(data_toplot, distance_metric));
    end
    
    if strcmp(method, 'tsne')
        fprintf('plot %d: %dms\n', t, time_ms(t))
        if t==1 || ~link_plots
            points = tsne(data_toplot,'Distance',distance_metric, 'NumPCAComponents', n_pc);
        else
            points = tsne(data_toplot,'Distance',distance_metric, 'InitialY', points, 'NumPCAComponents', n_pc);
        end
    elseif strcmp(method, 'mds')
        if t==1 || ~link_plots 
            points = mdscale(data_toplot, 2);
        else
            points = mdscale(data_toplot, 2, 'Start', points);
        end
    elseif strcmp(method, 'cmds')
        points = cmdscale(data_toplot, 2);
    elseif strcmp(method, 'pca')
        if link_plots
            tind = find(time_inds==t0);
            points = scores(((tind-1)*n_stim+1) : tind*n_stim,:);
        else
            [~,points] = pca(data_toplot, 'NumComponents', 2);
        end
    end
    
    if categories
        gscatter(points(:,1),points(:,2),cat_names(color_code(:,1)),[],[],[],'off');
    else
        scatter(points(:,1),points(:,2),[],color_code(:,1),'filled');
    end
    if size(color_code,2)==2
        hold on; scatter(points(:,1),points(:,2),45,color_code(:,2),'Linewidth',1.5); hold off
    end
    title(sprintf('%d ms',time_ms(t)))
end

colormap(jet);
if categories
    legend(cat_names,'Position',[0,0.6,0.1,0.1])
end
if link_plots
    linkaxes(ha,'xy');
end

end