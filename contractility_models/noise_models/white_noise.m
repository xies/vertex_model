function C = white_noise(ct,~,params)
%WHITE_NOISE
% Adds a white noise of given magnitude and bias.
% 
% USAGE: C = white_noise(ct,t,params)
%
% INPUT: ct - cell centroids (used to judge size)
%        t - time (not used)
%        params(1) - noise amplitude
%        params(2) - noise bias

C = params(1) * randn(1,size(ct,1)) + params(2);

end