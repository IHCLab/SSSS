function plot_result(d,Xhat_im)
%% B2 B3 B4 (10m) (true color)
load Yim_cell;
RGB= cat(3,Yim_cell{4},Yim_cell{3},Yim_cell{2});
fig_input_10m_B4B3B2= (RGB-min(min(min(RGB))))/(max(max(max(RGB)))-min(min(min(RGB))));
%% Input of SSSS
for i=1:12,
    Yim(:,:,i)= kron(Yim_cell{i},ones(d(i),d(i)));
    Yim(:,:,i)= (Yim(:,:,i)-min(min(Yim(:,:,i))))/(max(max(Yim(:,:,i)))-min(min(Yim(:,:,i))));
    Xhat_im(:,:,i)= (Xhat_im(:,:,i)-min(min(Xhat_im(:,:,i))))/(max(max(Xhat_im(:,:,i)))-min(min(Xhat_im(:,:,i))));
end
fig_input_60m_B1B9= cat(3,Yim(:,:,[1 10]),zeros(size(Yim(:,:,1)))); % B1 B9 (60m)
fig_input_20m_B5B6B7= Yim(:,:,[5 6 7]); % B5 B6 B7 (20m)
fig_input_20m_B8aB11B12= Yim(:,:,[9 11 12]); % B8a B11 B12 (20m)
%% Output of SSSS (super-resolved image)
fig1= Xhat_im(:,:,[1 10]); fig1= cat(3,fig1,zeros(size(fig1(:,:,1)))); % B1 B9 (60m)
fig2= Xhat_im(:,:,[5 6 7]); % B5 B6 B7 (20m)
fig3= Xhat_im(:,:,[9 11 12]); % B8a B11 B12 (20m)
%% plot
figure
subplot(1,7,1)
imshow(fig_input_10m_B4B3B2,'InitialMagnification','fit');
title('(a) RGB (10m)')
subplot(1,7,2)
imshow(fig_input_60m_B1B9,'InitialMagnification','fit');
title('(b) B1-9 (60m)')
subplot(1,7,3)
imshow(fig1,'InitialMagnification','fit');
title('(c) B1-9 (10m)')
subplot(1,7,4)
imshow(fig_input_20m_B5B6B7,'InitialMagnification','fit');
title('(d) B5-6-7 (20m)')
subplot(1,7,5)
imshow(fig2,'InitialMagnification','fit');
title('(e) B5-6-7 (10m)')
subplot(1,7,6)
imshow(fig_input_20m_B8aB11B12,'InitialMagnification','fit');
title('(f) B8a-11-12 (20m)')
subplot(1,7,7)
imshow(fig3,'InitialMagnification','fit');
title('(g) B8a-11-12 (10m)')
step= 0.01; % set position of the subplot
for i=2:7,
    im_now= subplot(1,7,i);
    P= get(im_now,'position');
    P(1)= P(1)-step*(i-1);
    set(im_now,'position',P)
end