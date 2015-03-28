function cells = create_squares(HEX_NUM_X, HEX_NUM_Y)
% Based on the parameters provided in we draw something in fourier space
% and then invert to get hexagon in normal space
%
%The output is an image matrix of hexagons


% arbitrary but ok
BIG_DIMENSION = 512; %px
BOUNDARY = 20; %px

if HEX_NUM_X > HEX_NUM_Y
    extx = BIG_DIMENSION;
    exty = round_n(extx * HEX_NUM_Y / HEX_NUM_X, 2);
else
    exty = BIG_DIMENSION;
    extx = round_n(exty * HEX_NUM_X / HEX_NUM_Y, 2);
end
ax = round_n(extx - BOUNDARY, 2);
ay = round_n(exty - BOUNDARY, 2);
    
% extx = 512;
% exty = 256;
% ax = 480;
% ay = 210;

% the spatial frequency -> number of cells
% this is arbitrary but at least resembles reasonable numbers
rr = max(HEX_NUM_X, HEX_NUM_Y);

% size of the image to be fft'd
imgsz = max(extx, exty);

% rotate by pi/2

% create the fourier transform of hexagons
q = zeros(imgsz);
for i=0:3
    q(round(imgsz/2 + 1 + rr*sin(i*pi/2 + pi/4)), ...
        round(imgsz/2 + 1 + rr*cos(i*pi/2 + pi/4))) = 1;
end
%

% inverse fourier transform
th = 0.5;
qq = -abs(ifft2(fftshift(q)));

qq = qq-min(qq(:));
qq = qq/max(qq(:));
qq = im2bw(qq-th,0);
qq = bwmorph(qq,'skel',Inf);

cells = zeros(exty,extx);
cells(exty/2-ay/2+1:exty/2+ay/2,extx/2-ax/2+1:extx/2+ax/2) = ...
   qq(imgsz/2-ay/2+1:imgsz/2+ay/2,imgsz/2-ax/2+1:imgsz/2+ax/2);

end
