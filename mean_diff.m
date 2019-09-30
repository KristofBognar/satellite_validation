function [ mean_diff, sigma, sigma_mean ] = mean_diff( type, x1, x2, rmsd )
%[mean_diff,sigma,sigma_mean]=mean_diff(type, x1, x2 ) Calculates mean absolute
%or relative difference
%
% INPUT:    type: 'abs' or 'rel', for absolute and relative difference. Relative 
%                   difference is taken w.r.t the mean of the two datapoints
%                 'rmsd_abs' or 'rmsd_rel', for root-mean-square deviation
%                   using absolute or relative differences as input
%           x1, x2: coincident measurements (must be same size)
%
% OUTPUT:   mean_diff: mean absolute or relative difference (rel. diff. in percent)
%           sigma: standard deviation of differences
%           sigma_mean: error on the mean =sigma/sqrt(N)
%   
% Kristof Bognar, December 2017

%% consistency checks
if length(x1)~=length(x2), error('Input arrays must be the same size'); end

% if any(isnan(x1+x2)), error('Input data contains NaNs'); end

%% calculate values

% get differences between individual pairs (absolute or relative)
[ diff ] = rel_abs_diff( type, x1, x2 );

if nargin==3
    % if only mean differences are needed
    mean_diff=nanmean(diff);

    % calculate std and standard error
    sigma=std(diff);
    sigma_mean=sigma/sqrt(length(x1));
    
elseif strcmp(rmsd,'rmsd')
    % RMSD values are needed
    mean_diff=sqrt(nanmean(diff.^2));

    % no errors returned
    sigma=NaN;
    sigma_mean=NaN;

end

end

