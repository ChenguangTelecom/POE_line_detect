close all
clear all
clc

%type mex poe_line_detect.c poe.c

BW=zeros(512);
BW(256,:)=1.0;%BW should be the binary edge map of the image
BW=imread('P1040833.png');
BW=double(BW)/255.0; % binary edge maps and with the value of 0 and 1
tic
W=10.0; %radius of the POE kernel
P=16.0; % number of orientations in POE
tau=pi/16; % angle tolerance in the region growing of POE
lines=poe_line_detect(BW,W,P,tau);

size(lines)
toc
figure,imshow(BW)
figure,imshow(zeros(size(BW)))
hold on
for k=1:size(lines,1)
xy=[lines(k,2),lines(k,1);lines(k,4),lines(k,3)];
plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
end
