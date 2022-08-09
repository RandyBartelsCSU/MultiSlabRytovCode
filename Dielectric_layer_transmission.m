clc; close all; clear

lambda      = 1;                  % central wavelength (microns)
n_imm       = 1; %1.367;                % refractive index of immersion media
n_o         = 2;                 % RI of the object
k0=(2*pi)/lambda;
k=k0*n_imm;

ps          = .05;                 % pixel size (x,y,z) in object space (microns)
N           = 2^6;                  % lateral pixel dimension 
x           = ps*(-N/2:N/2-1);      % 1D axis in x
y=x;
dfx         = 1/(N*ps);             % Fourier spacing of padded axis
[X, Y]=meshgrid(x,y);
fx          = dfx*(-N/2:N/2-1);     % 1D axis in fx
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

fxx         = ifftshift(fxx);       % FFT shifting Fourier axes
fyy         = ifftshift(fyy);       % FFT shifting Fourier axes

L = 10;
Nz=L/ps;
z=ps*(0.5:1:Nz-0.5);

D = 3;
Noz=D/ps;

RI = n_imm*ones(N,N,Nz);
RI(:,:,1:Noz)=n_o;

V=-(k0)^2*((RI).^2-n_imm^2);
eps=0;

U_in = exp(-1i*k*0.5*ps)*ones(N,N);

EB=MultiLayerBorn(fxx,fyy,lambda,n_imm,ps,V,U_in,eps,'Vol');%.017
EB=EB(:,:,2:end);

ER=MultiLayerRytov(fxx,fyy,lambda,n_imm,ps,V,U_in,eps,'Vol');
ER=ER(:,:,2:end);

figure
subplot(4,2,1)
imagesc(z,x,squeeze(V(:,32,:)))
colorbar
title('V')
subplot(4,2,3)
imagesc(z,x,squeeze(angle(EB(:,32,:))))
colorbar
title('\angle E_{MLB}')
subplot(4,2,4)
imagesc(z,x,squeeze(abs(EB(:,32,:))))
colorbar
title('|E_{MLB}|')
subplot(4,2,5)
imagesc(z,x,squeeze(angle(ER(:,32,:))))
colorbar
title('\angle E_{MSR}')
subplot(4,2,6)
imagesc(z,x,squeeze(abs(ER(:,32,:))))
colorbar
title('|E_{MSR}|')
subplot(4,2,7)
plot(z,squeeze(angle(EB(32,32,:))))
colorbar
hold on
plot(z,squeeze(angle(ER(32,32,:))))
legend('MLB','MSR')
title('\angle E')
subplot(4,2,8)
plot(z,squeeze(abs(EB(32,32,:))))
hold on
plot(z,squeeze(abs(ER(32,32,:))))
legend('MLB','MSR')
title('|E|')

