function fig=rsa_plot_rdm(rdm, time_ms, varargin)
% Plot representational dissimilarity matrices (RDM) corresponding to the
% given time-points (given in ms).
%
% Default settings for conversion of ms to indices in RDM:
%   sampling_rate   = 250 Hz
%   baseline        = 100 ms
% (to change see optional input)
%
% Example usage:
%   - fig = rsa_plot_rdm(rdm, [-50:50:300])
%       Plot the RDMs corresponding to 50ms before stimulus onset to 300ms
%       after onset in jumps of 50ms (8 plots). Conversion to indices is
%       done using default settings.
%   - fig = rsa_plot_rdm(rdm, [-50:50:300], 'categories', cat_inds, cat_names)
%       Order the stimuli according to cat_inds, and mark category
%       boundaries in the rdm.
%   - fig = rsa_plot_rdm(rdm, [], 'sampling_rate', 512)
%       Plot RDMs at the default time-points (see below), and assume
%       sampling rate is 512 Hz for conversion to indices.
% 
% Input:
%   rdm (RDMs per time-point, dimensions: n_stim * n_stim * seg_duration)
%   time_ms (Optional): plot the RDMs corresponding to these time-points
%       (in ms around stimulus onset, e.g. 50ms before onset is -50, and
%       50ms after stimulus onset is +50).
%       Default conversion to indices in rdm is done using sampling rate
%       250 Hz and 100ms baseline. Change using optional arguments.
%       If no time_ms is supplied the default is 6 evently spaced time-points
%       where the distance from segment boundaries is 1/2 distance between
%       time-points, e.g. if the segment is 300 samples the default is to
%       plot samples 25, 75, 125, 175, 225 and 275.
%   Additional optional arguments (varargin):
%       - To order the RDM according to categories use: 'categories'
%           followed by (in this order!):
%           1) cat_inds - a numeric vector (length n_stim) indicating
%              the category assignment of each stimulus. Use numbers from 1
%              to n_categories.
%           2) cat_names - a cell array of the name of each category.
%       - 'sampling_rate' followed by a number (in Hz) - change the
%           sampling rate for the conversion of ms to indices.
%       - 'baseline' followed by a number (in ms) - change the baseline
%           duration for the conversion of ms to indices.
%       - 'cross_validated' if the distances are cross validated and you
%           don't want to remove the diagonal
%
% Output: figure with RDMs. Distances between subplots and location of
% category numbers might require tweeking for your computer's settings.
%
% Written by Gal Vishne, lab of Leon Y. Deouell, 2019-2020
% Send bug reports and requests to gal.vishne@gmail.com

n_stim = size(rdm,1); duration = size(rdm,3);

categories = false;
baseline = 100;         % in ms
sampling_rate = 250;    % in Hz
cross_validated = false;
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
        case 'cross_validated'
            cross_validated = true;
            arg = arg + 1;
        otherwise
            error(['Unknown optional argument name: ' varargin{arg} '.']);
    end
end

if exist('time_ms','var')
    if isempty(time_ms)
        time_inds = linspace(1, duration, 13);
        time_inds = round(time_inds(2:2:13));
        time_ms = ((time_inds-1)*1000/sampling_rate - baseline);
    else
        time_inds = floor((time_ms+baseline)*sampling_rate/1000)+1; %turn ts into indices
        if min(time_inds) < 1 || max(time_inds) > duration
            error('Given temporal range exceeds segment duration');
        end
    end
else
    time_inds = linspace(0, duration, 13);
    time_inds = round(time_inds(2:2:13));
    time_ms = ((time_inds-1)*1000/sampling_rate - baseline);
end

% don't include the diagonal in the color range
all_dists_plotted = rdm(:,:,time_inds);
if ~cross_validated
    for t = 1:length(time_inds)
        all_dists_plotted(:,:,t)=all_dists_plotted(:,:,t)+diag(nan(size(all_dists_plotted,1),1));
    end
end
color_range = [nanmin(all_dists_plotted(:)), nanmax(all_dists_plotted(:))];

if length(time_inds)<3
    num_rows=1;
else
    num_rows=2;
end
num_cols=round(length(time_inds)/num_rows);
fig = figure('Units','Normalized','Position',[0.3 0.3 0.2+0.08*num_cols 0.4]);

if categories
    [cat_inds_sorted, sorting_order] = sort(cat_inds);
    for t = 1:length(time_inds)
        all_dists_plotted(:,:,t) = all_dists_plotted(sorting_order, sorting_order, t);
    end
    tick_locations = find([boolean(0);diff(cat_inds_sorted)~=0])-0.5;
    categ_locations = [tick_locations(1)/2;tick_locations+[diff(tick_locations);n_stim-tick_locations(end)]/2];
    categ_nums = num2str((1:length(cat_names))');
end

for t=1:length(time_inds)
    w_inset = [0.12 0.02]; h_inset = [0.03 0.07]; w_gap = 0.018; h_gap = 0.08;
    width = (1-sum(w_inset)-(num_cols-1)*w_gap)/num_cols; height = (1-sum(h_inset)-h_gap)/num_rows;
    if t<=num_cols
        pos = [w_inset(1)+(t-1)*(width + w_gap) h_inset(1)+(h_gap+height)*(num_rows-1) width height];
    else
        pos = [w_inset(1)+(t-1-num_cols)*(width + w_gap) h_inset(1) width height];
    end
    axes('Position',pos)

    rdm_toplot = all_dists_plotted(:,:,t);
    imagesc(rdm_toplot,color_range);hold on
    tit=title(sprintf('%d ms',time_ms(t)));
    if categories
        txt1 = text(categ_locations,zeros(size(categ_locations))-6,categ_nums,'HorizontalAlignment','center');
        txt2 = text(zeros(size(categ_locations))-4,categ_locations,categ_nums,'HorizontalAlignment','center');
        for i=tick_locations'
            plot([0,n_stim],[i,i],'k');
            plot([i,i],[0,n_stim],'k');
        end
        set(tit,'position',get(tit,'position')-[0 6 0])
        set(gca,'XAxisLocation','top','XTick',tick_locations,'XTickLabels','',...
            'YTick',tick_locations,'YTickLabels','');
    else
        set(gca,'XAxisLocation','top','XTickLabels','','YTickLabels','');
    end
    colormap jet
end
cb=colorbar;
set(cb,'position',[.04 .1 .03 .3])

if categories
    txt='Categories:';
    for i=1:length(categ_nums)
       txt = [txt newline categ_nums(i) '-' cat_names{i}];
    end
    annotation('textbox',[0.01,0.8,0.1,0.1],'String',txt,'FitBoxToText','on','Margin',4);
end
hold off

end