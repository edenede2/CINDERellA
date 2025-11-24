function [pa_loc,postPend1,n1_ls]=drawpa3(pb_H0,n1_Spa_kj1_kj2,postPend1,u)
% DRAWPA3 - Search move destination specified by random number u
% This function is used during MCMC sampling to select parent sets
%
% Inputs:
%   pb_H0         - Probability vector
%   n1_Spa_kj1_kj2 - Local scores for parent sets
%   postPend1     - Posterior pending value
%   u             - Random number for selection
%
% Outputs:
%   pa_loc        - Selected parent location index
%   postPend1     - Updated posterior pending value
%   n1_ls         - Local score for selected parent set

n1_ls=[];
pa_loc=[];
postP=postPend1+cumsum(pb_H0);

if isempty(postP)
    return;
end

if postP(end)<u
    postPend1=postP(end);
else
    if postP(1)>u
        pa_loc = 1;
    else
        pa_loc = find(postP <u,1,'last')+1;
        postPend1=postP(pa_loc-1);
    end
    
    n1_ls = n1_Spa_kj1_kj2(pa_loc);
end
