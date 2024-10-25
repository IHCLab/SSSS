%=====================================================================
% Programmer: Chia-Hsiang Lin (Steven)
% Web: https://sites.google.com/view/chiahsianglin/
% E-mail: chiahsiang.steven.lin@gmail.com
% Date: Nov. 13, 2019
%======================================================================
% Reference:
% C.-H. Lin, and J. M. Bioucas-Dias,
% ``An Explicit and Scene-adapted Definition of Convex Self-similarity Prior with Application to Unsupervised Sentinel-2 Super-resolution,"
% IEEE Transactions on Geoscience and Remote Sensing, 2019.
%======================================================================
% Sentinel-2 Superresolution via Scene-adapted Self-similarity (SSSS)
% [Xim,time] = SSSS(Yim,rv,dx,dy,sdf,lambda,mu)
%----------------------------------------------------------------------
% INPUT
% Yim: nr*nc*nb multi-resolution data with 0 inserted
% rv: nb*1 vector with entries in {1,2,6}
% (dx,dy,sdf) define the blurring kenel 'fspecial('gaussian',[dx,dy],sdf(i))' of band i
% lambda: regularization parameter
% mu: ADMM penalty parameter
%----------------------------------------------------------------------
% OUTPUT
% Xim: superresolved image whose dimension is the same as Yim
% time: computational time
%----------------------------------------------------------------------
% Assumption
% 1. dx and dy are odd numbers (so that kenel has a center)
% 2. nr and nc can be even or odd, but should be multiples of 6 (hence even)
% 3. blurring matrix is BCCB
% 4. sampling at the top-left pixel
%========================================================================
function [Xim,time]= SSSS(Yim,rv,dx,dy,sdf,lambda,mu)
%% setting
t0= clock;
degreeNET= 1; % degree of network (i.e., #(edges) per patch)
patch_size= 6;
band_for_net_learn= [2 3 4 8]; % 2 (blue), 3 (green), 4 (red), or 8 (NIR)
[nr,nc,nb]= size(Yim);
n= nr*nc;
boundary_width= 6; % (dx-1)/2 = (dy-1)/2 = 6
Y= (reshape(Yim,n,nb))';
p= 5; % p = nb
SubspaceLearn= 1; % p < nb
if (SubspaceLearn==1), p= 5; end
epoch= 30; % number of iterations between 2 subspace learning
num_ADMM_itr= 100;
%% blurring kernel & blurring matrix
middlel= 1+round(nr/2); % dx and dy are assumed odd
middlec= 1+round(nc/2); % dx and dy are assumed odd
for i=1:nb,
    if sdf(i) > 0,
    K= fspecial('gaussian',[dx,dy],sdf(i)); % cropped kernel (dx*dy)
    KK= zeros(nr,nc);
    KK(middlel-(dy-1)/2:middlel+(dy-1)/2,middlec-(dx-1)/2:middlec+(dx-1)/2)= K;
    KK= fftshift(KK);
    KK= KK/sum(sum(KK));
    KK(:,:,i)= KK; % full kernel (nr*nc), such that K centered at (1,1)
    else
    KK(:,:,i)= zeros(nr,nc);
    KK(1,1,i)= 1; % full kernel (nr*nc), such that K centered at (1,1)
    end
    FBM(:,:,i)= fft2(KK(:,:,i));
end
%% network learning
[Edge,Alpha]= NetLearn(mean(Yim(:,:,band_for_net_learn),3),patch_size,boundary_width,degreeNET);
NoE= length(Alpha); % number of edges
%% operators
ConvCM = @(X,FKM) reshape(real(ifft2(fft2(reshape(X',nr,nc,size(X,1))).*FKM)),nr*nc,size(X,1))';
DH1 = @(X) reshape(sum(im2col(X,[nr/1,nc/1],'distinct'),2),nr/1,nc/1); % sum the r^2 subimages in X (r=1)
DH2 = @(X) reshape(sum(im2col(X,[nr/2,nc/2],'distinct'),2),nr/2,nc/2); % sum the r^2 subimages in X (r=2)
DH6 = @(X) reshape(sum(im2col(X,[nr/6,nc/6],'distinct'),2),nr/6,nc/6); % sum the r^2 subimages in X (r=6)
D1 = @(X) repmat(X,1,1); % copy an image r^2 times (r=1)
D2 = @(X) repmat(X,2,2); % copy an image r^2 times (r=2)
D6 = @(X) repmat(X,6,6); % copy an image r^2 times (r=6)
%% ADMM
d= zeros(nb*n,1);
v= d;
Dij= zeros(n,p,NoE);
Vij= Dij;
for i=1:nb,
    Xim(:,:,i)= PAN4(Yim(:,:,i),Yim(:,:,[2 3 4 8]),rv(i),FBM(:,:,i));
end
UZ= (reshape(Xim,n,nb))'; % UZ=U*Z
[U,~,~]= svd(UZ*UZ');
if p < nb,
    U= U(:,1:p);
    Z= U'*UZ;
else
    U= eye(nb);
    Z= UZ;
end
UTU_plus_I_inv= inv( (U')*U+NoE*eye(p) );
FBMC= conj(FBM);
BTY= ConvCM(Y,FBMC);
for b=1:nb,
    S= FBM(:,:,b); % diagnal of \Sigma, in image format
    r= rv(b);
    if (r==1), IF{b}= 1./(r^2+(DH1(abs(S.*S)))/mu); end
    if (r==2), IF{b}= 1./(r^2+(DH2(abs(S.*S)))/mu); end
    if (r==6), IF{b}= 1./(r^2+(DH6(abs(S.*S)))/mu); end
end
for ij=1:NoE,
    [idx_pi_pj{ij},PTP_plus_I_inv{ij}]= patch_idx(nr,nc,Edge(ij,:),patch_size,mu,lambda,Alpha(ij)); % for Vij update
end
for i=1:num_ADMM_itr,
    NU= U*Z-(reshape(d,n,nb))';
    AUX= BTY+mu*NU;
    for b=1:nb,
        input= reshape(AUX(b,:),[],1);
        r= rv(b);
        S= FBM(:,:,b); % diagnal of \Sigma, in image format
        if (r==1), vv= ifft2(D1(DH1(fft2(reshape(input,nr,nc)).*S).*IF{b}).*conj(S)); end
        if (r==2), vv= ifft2(D2(DH2(fft2(reshape(input,nr,nc)).*S).*IF{b}).*conj(S)); end
        if (r==6), vv= ifft2(D6(DH6(fft2(reshape(input,nr,nc)).*S).*IF{b}).*conj(S)); end
        v(1+(b-1)*n:b*n)= input/mu-real(reshape(vv,[],1))/mu/mu;
    end
    % start Vij update
    for ij=1:NoE,
        Delta= Z'-Dij(:,:,ij);
        Vij(:,:,ij)= Delta; % entries not in the two patches
        Vij(idx_pi_pj{ij},:,ij)= PTP_plus_I_inv{ij}*((mu/lambda/Alpha(ij))*Delta(idx_pi_pj{ij},:)); % entries in the two patches
    end
    % start z update
    Z= ((((reshape(v+d,n,nb))*U)+(sum(Vij+Dij,3)))*UTU_plus_I_inv)';
    % start dual update
    d= d+v-reshape((U*Z)',[],1);
    Dij= Dij+Vij-repmat(Z',[1 1 NoE]);
    % start subspace learning
    if ((SubspaceLearn==1)&(i==epoch)),
        term2= Z*BTY';
        for b=1:nb,
            r= rv(b);
            BZT= ConvCM(Z,repmat(FBM(:,:,b),[1 1 p]));
            BZT3D= reshape(BZT',[nr nc p]);
            MBZT3D= zeros(size(BZT3D));
            MBZT3D(1:r:end,1:r:end,:)= BZT3D(1:r:end,1:r:end,:);
            MBZT= sparse(reshape(MBZT3D,[n p]));
            term1= inv(MBZT'*MBZT);
            U(b,:)= (term1*term2(:,b))';
        end
        UTU_plus_I_inv= inv((U')*U+NoE*eye(p));
    end
end
Xim= reshape(v,[nr nc nb]);
time= etime(clock,t0);
%% subprogram
function SYN= PAN4(Y,X4,r,FBMi)
if (r==0),
    SYN= Y;
else
    [nr,nc]= size(Y);
    n= nr*nc;
    ConvCM = @(X,FKM) real(ifft2(fft2(reshape(X',nr,nc,size(X,1))).*FKM));
    BX4= ConvCM((reshape(X4,n,4))',repmat(FBMi,[1 1 4]));
    d= reshape(Y(1:r:end,1:r:end),[],1);
    C= reshape(BX4(1:r:end,1:r:end,:),[],4);
    coef4= lsqnonneg(C,d);
    SYN= reshape((reshape(X4,[],4))*coef4,[nr nc]);
end
return;
%% subprogram
function [idx_pi_pj,PTP_plus_I_inv]= patch_idx(nr,nc,edge,ps,mu,lambda,alpha)
i= edge(1);
j= edge(2);
n= nr*nc;
yi= mod(i,nr); if (yi==0), yi=nr; end
xi= 1+((i-yi)/nr);
yj= mod(j,nr); if (yj==0), yj=nr; end
xj= 1+((j-yj)/nr);
DH = @(X) reshape(sum(im2col(X,[nr,nc],'distinct'),2),nr,nc); % sum the r^2=4 subimages in X
tempi= zeros(2*nr,2*nc);
tempi(yi:yi+ps-1,xi:xi+ps-1)= ones(ps,ps);
mapi= DH(tempi);
idx_pi= find(mapi==1);
Pi= sparse(zeros(ps*ps,n));
Pi(:,idx_pi)= speye(ps*ps);
tempj= zeros(2*nr,2*nc);
tempj(yj:yj+ps-1,xj:xj+ps-1)= ones(ps,ps);
mapj= DH(tempj);
idx_pj= find(mapj==1);
Pj= sparse(zeros(ps*ps,n));
Pj(:,idx_pj)= speye(ps*ps);
idx_pi_pj= find((or(mapi,mapj))==1);
cardI= length(idx_pi_pj);
Pi_minus_Pj= sparse(Pi-Pj);
Pi_minus_Pj_active= sparse(Pi_minus_Pj(:,idx_pi_pj));
PTP_plus_I_inv= inv(((Pi_minus_Pj_active')*Pi_minus_Pj_active)+((mu/lambda/alpha)*speye(cardI)));
return;
%% subprogram
function [Edge,Alpha]= NetLearn(IM,ps,bw,degree)
[nr,nc]= size(IM);
% find overlapping patches
temp= repmat(IM,2,2); % circulant
patches= im2col(temp(1:nr+ps-1,1:nc+ps-1),[ps,ps],'sliding'); % patches(1,:) == reshape(IM,1,[])
% find invalid patch index (at the image boundary)
all_idx= ones(nr,nc);
all_idx(1+bw:nr-bw-(ps-1),1+bw:nc-bw-(ps-1))= 0; % mark valid (interior) patches by 0
invalid_idx= find(all_idx==1);
% for each non-overlapping patch, find its most similar patch
k=0;
for i= 1+bw:ps:nr-bw-(ps-1),
    for j= 1+bw:ps:nc-bw-(ps-1),
        k= k+1;
        idx= i+(j-1)*nr;
        diff= sum(((patches(:,idx)*ones(1,nr*nc))-patches).^2,1);
        diff(idx)= inf;
        diff(invalid_idx)= inf; % do not consider boundary patches
        [dist(k,1),idx_similar]= min(diff);
        Edge(k,:)= [idx,idx_similar];
        if (degree>1),
            for edge_idx=2:degree,
                k= k+1;
                diff(idx_similar)= inf;
                [dist(k,1),idx_similar]= min(diff);
                Edge(k,:)= [idx,idx_similar];
            end
        end
    end
end
dist= sqrt(dist);
Alpha= dist.^(-1);
return;