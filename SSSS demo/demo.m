clear; close all; clc
%% upsampling matrix
load Yim_cell; % real Sentinel-2 data
rv= [6 1 1 1 2 2 2 1 2 6 2 2]';
for i=1:12,
    topleft= zeros(rv(i),rv(i)); topleft(1,1)= 1;
    Yim(:,:,i)= kron(Yim_cell{i},topleft);
end
[C,D,nb]= size(Yim);
%% blurring matrix
dx= 13; dy= 13; % kernel filter support
limsub= 6; % due to BCCB convolution, one needs remove border (dx-1)/2
mtf= [.32 .26 .28 .24 .38 .34 .34 .26 .23 .33 .26 .22];
sdf= rv.*sqrt(-2*log(mtf)/pi^2)';
sdf(rv==1)= 0; % (dx,dy,sdf) define the blurring kenel 'fspecial('gaussian',[dx,dy],sdf(i))' of band i
%% algorithm
t0= clock;
lambda= 0.1; % regularization parameter in SSSS
mu= 0.1; % ADMM penalty parameter in SSSS
nr= 72; A= (C-2*limsub)/(nr-2*limsub);
% nc= 108; B= (D-2*limsub)/(nc-2*limsub);
row_start= 1; row_end= row_start+nr-1
[Xhat_im(1:nr,:,:),time]= SSSS(Yim(row_start:row_end,:,:),rv,dx,dy,sdf,lambda,mu);
for i=2:A,
    row_start= row_start+nr-2*limsub; row_end= row_start+nr-1
    [temp,time]= SSSS(Yim(row_start:row_end,:,:),rv,dx,dy,sdf,lambda,mu);
    Xhat_im(row_start+limsub:row_end,:,:)= temp(1+limsub:end,:,:);
end
% Xhat_im1= Xhat_im(limsub+1:end-limsub,limsub+1:end-limsub,:); % merged & border removed
TIME= etime(clock,t0)
%% plot
plot_result(rv,Xhat_im)