function vid = movie(tissues,varargin)
% Make a movie of tissue evolving. Can return just the
% outlines of cells, or also shade-in the active cells.
%
% USAGE: I = movie(tisues);
%        I = movie(tisues,'showActive');

num_frames = numel(tissues);
if nargin == 1, opt = {'none'};
else
    opt = varargin;
end

ind = find( strcmpi(varargin,'save') );
if ~isempty(ind)
    OUT_PATH = varargin{ind+1};
else
    OUT_PATH = '~/Desktop/model.avi';
end

vid = VideoWriter( OUT_PATH );
vid.Quality = 100;
vid.FrameRate = 7;
open(vid);
for f = 1:num_frames
    tissues(f).draw(opt{:});
    title(['Time = ' num2str(tissues(f).t)]);
    vid.writeVideo(getframe(gcf));
end
close(vid);
end