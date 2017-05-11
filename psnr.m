function psnr = psnr(x,xr);
%    psnr = psnr(x,xr);
%    Compute the Peak signal-noise ratio (PSNR) of the 
%    reconstructed signal xr comparing to the original
%    signal x

xmax=max(max(x));
% For 8-bit grayscale image, xmax is usually set to 255!

error=abs(x-xr);
mse=mean(mean(error.*error));
psnr=10*log10(xmax*xmax/mse);
disp(['Mean square error (MSE)       : ',num2str(mse)]);
disp(['Peak signal-noise ratio (PSNR): ',num2str(psnr),' dB']);
