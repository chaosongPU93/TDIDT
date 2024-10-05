% img2video 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will read in several Bitmap type image files into frames, and
% convert them into videos
% 
% --There is an instability in Linux version matlab. When the iterative plots
%   are first made, the fig frame in the first iteration is not the exactly 
%   same as the rest ones, even though your code is the same, such that the 
%   saved figures have a slightly different size. The different sizes lead
%   to errors when converting into videos which require exactly the same frame
%   size, no matter on which system are you making the video. To get around
%   it, run the joint inversion again with the iterative plot window open 
%   so that the figure size will remain constant. 
% --Linux version matlab does NOT support video types other than 'avi'. But
%   it will distort the shape of images, no matter what. Sometimes it will be
%   normal, but i have no idea about how, so use MacOS version to run this code
%   whenever possbile
% --MacOS version supports other types, including 'mp4'. MP4 format is just 
%   better than AVI, more compact in size.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/02/08
% Last modified date:   2022/02/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
% clear
if exist('vid','var')
  clear vid
end

%Image directory
imgdir = '/home/chaosong/Pictures/';
imgfiles = dir(fullfile(imgdir,'IndepIter*.jpg'));
% imgfiles = dir(fullfile(imgdir,'JointIter*.jpg'));

%Change video file name and file type according to the operating system
vidfile = 'indep.avi';
% vidfile = 'joint.avi';
% vidfile = 'itersteps.mp4';

vidprofile = 'Motion JPEG AVI';
% vidprofile = 'MPEG-4';

vid = VideoWriter(fullfile(imgdir,vidfile),vidprofile);
vid.FrameRate = 1;  
vid.Quality = 100;

open(vid);
for i = 1: size(imgfiles,1)
  fname = fullfile(imgdir,strcat('IndepIter',num2str(i),'.jpg'));
%   fname = fullfile(imgdir,strcat('JointIter',num2str(i),'.jpg'));
  img = imread(fname);  % read image
  frame = im2frame(img);  % image to frame
  writeVideo(vid, frame); % write frame into the video
%   writeVideo(vid, img); % write frame into the video
end        
close(vid);








