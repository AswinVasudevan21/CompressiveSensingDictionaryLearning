clear all
clc
close all

addpath('./Utils/rwt-master/src');
addpath('./Data');
addpath('./Recon');

% load image
img = imread('boat256.tif');
img = double(rgb2gray(img));
img = img./max(img(:));
[DIM1,DIM2] = size(img);

n = 8; % patch size
N = n^2; % vectorization later
K = 10; % sparsity
lambda = 0.1; % regularization parameter for tradeoff btw data fitting and l1 norm
M = round(3*K); % number of measurments
h = daubcqf(4,'min'); % setting for wavelet filter
L = log2(N); % setting for wavelet filter (level)

param.L = K;
noise = 0;
img_recon = zeros(size(img));
for i = 1:n:DIM1
    for j = 1:n:DIM2
        patch = img(i:i+n-1,j:j+n-1);
        patch = patch(:);
        
        %% sensing matrix Phi
        % Gaussian
        %Phi = randn(M,N);
        %matrixNorm = Phi.'*Phi;
        %matrixNorm = sqrt(diag(matrixNorm)).';
        %Phi = Phi./repmat(matrixNorm, [M,1]);
        % subsampling matrix
        D = dctmtx(N);
        iD = D^-1;
        perm = randperm(N);
        Phi = iD(perm(1:M),:);
        %% measurement vector
        y = Phi*dct(patch); %dct
%       y = Phi*mdwt(patch,h,L); %dwt
        
        y = y + randn(size(y))*noise;
        
        %% test various reconstruction methods here
%         alpha = OMP(y, Phi, K); % OMP
%         alpha = SP(y,Phi,K); % SP
%         alpha = FISTA(Phi,y); % FISTA
        alpha = ALM(Phi,y); % ALM
%       [alpha,~,~,~,~,~] = GPSR_Basic(y,Phi,lambda); %GPSR

        patch_recon = reshape(idct(alpha),n,n); % idct
%       patch_recon = reshape(midwt(alpha,h,L),n,n); % idwt
        img_recon(i:i+n-1,j:j+n-1) = patch_recon;
    end
end;
figure;imshow(img,[]);
figure;imshow(img_recon,[]);
psnr1 = psnr(img,img_recon);
