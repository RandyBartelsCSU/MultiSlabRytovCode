function [RI,sdev,avg]=RIGenerator3D(tissuemodel,l,sigma,nu,x,y,z,precision)

% This function generates RI material based on several different
% correlation function models

%Inputs:
% tissuemodel: string; tissue model is either 'MW'(Mattle-Whittern), 'BG'(Booker-Gordon), or 'Gauss'(Gaussian) 
% l: correlation lengh in microns
% sigma: standard deviation of desired RI
% nu: parameter for MW...
% x,y,z spatial grids(vectors) in microns for RI
switch precision
    case 'single'
x=single(x);
y=single(y);
z=single(z);
Nx=max(size(x));
Ny=max(size(y));
Nz=max(size(z));
% ps=mean(diff(x));
% psz=mean(diff(z));
% 
% dfx=1/(Nx*ps);
% dfy=1/(Ny*ps);
% dfz=1/(Nz*psz);
% fx=dfx*[-Nx/2:Nx/2-1];
% fy=dfy*[-Ny/2:Ny/2-1];
% fz=dfz*[-Nz/2:Nz/2-1];
% [Fxx,Fyy,Fzz]= meshgrid(fx,fy,fz); 
% Fxx= ifftshift(Fxx);       % FFT shifting Fourier axes
% Fyy= ifftshift(Fyy);       % FFT shifting Fourier axes
% Fzz= ifftshift(Fzz);       % FFT shifting Fourier axes
%[X,Y,Z]=meshgrid(x,y,z);
R=zeros(length(x),length(y),length(z),'single');  % important to create R this way instead of meshgrid because meshgrid uses 3X the data (X Y and Z) 
for ii=1:length(x)
    for jj=1:length(y)
        for kk=1:length(z)
            R(ii,jj,kk)=sqrt(x(ii).^2+y(jj).^2+z(kk).^2);
        end
    end
end
%% Generate random field
Mean=0;
WN3=1*randn(Nx,Ny,Nz,'single')+Mean; % variance of 1 (also std)
%% Correlation Functions

% Mattle-Whittern:
BMW3=@(r) ((sigma.^2)./((2.^(nu-1)).*gamma(nu))).*((r./l).^nu).*besselk(nu,r./l);
% Booker-Gordon:
BBG3=@(r) (sigma.^2)*exp(-r/l);
% Gaussian:
BGauss3=@(r) (sigma.^2)*exp(-(r/l).^2);
%% Calculate RI
switch tissuemodel
    case 'MW'
        RI=ifftn(fftn(WN3).*sqrt(abs(fftn(BMW3(R)))));
    case 'BG'
        RI=ifftn(fftn(WN3).*sqrt(abs(fftn(BBG3(R)))));
    case 'Gauss'
        RI=ifftn(fftn(WN3).*sqrt(abs(fftn(BGauss3(R)))));
    otherwise
        disp('Model is either MW(Mattle-Whittern), BG(Booker-Gordon), or Gauss(Gaussian)')
end
%% Ensure stats match up
%V=-(k)^2*((RI+n_imm).^2-n_imm^2);
sdev=std(RI,[],'all');
avg=mean(RI,'all');
    case 'double'
        
Nx=max(size(x));
Ny=max(size(y));
Nz=max(size(z));
% ps=mean(diff(x));
% psz=mean(diff(z));
% 
% dfx=1/(Nx*ps);
% dfy=1/(Ny*ps);
% dfz=1/(Nz*psz);
% fx=dfx*[-Nx/2:Nx/2-1];
% fy=dfy*[-Ny/2:Ny/2-1];
% fz=dfz*[-Nz/2:Nz/2-1];
% [Fxx,Fyy,Fzz]= meshgrid(fx,fy,fz); 
% Fxx= ifftshift(Fxx);       % FFT shifting Fourier axes
% Fyy= ifftshift(Fyy);       % FFT shifting Fourier axes
% Fzz= ifftshift(Fzz);       % FFT shifting Fourier axes
%[X,Y,Z]=meshgrid(x,y,z);
R=zeros(length(x),length(y),length(z));
for ii=1:length(x)
    for jj=1:length(y)
        for kk=1:length(z)
            R(ii,jj,kk)=sqrt(x(ii).^2+y(jj).^2+z(kk).^2);
        end
    end
end
%% Generate random field
Mean=0;
WN3=1*randn(Nx,Ny,Nz)+Mean; % variance of 1 (also std)
%% Correlation Functions

% Mattle-Whittern:
BMW3=@(r) ((sigma.^2)./((2.^(nu-1)).*gamma(nu))).*((r./l).^nu).*besselk(nu,r./l);
% Booker-Gordon:
BBG3=@(r) (sigma.^2)*exp(-r/l);
% Gaussian:
BGauss3=@(r) (sigma.^2)*exp(-(r/l).^2);
%% Calculate RI
switch tissuemodel
    case 'MW'
        RI=ifftn(fftn(WN3).*sqrt(abs(fftn(BMW3(R)))));
    case 'BG'
        RI=ifftn(fftn(WN3).*sqrt(abs(fftn(BBG3(R)))));
    case 'Gauss'
        RI=ifftn(fftn(WN3).*sqrt(abs(fftn(BGauss3(R)))));
    otherwise
        disp('Model is either MW(Mattle-Whittern), BG(Booker-Gordon), or Gauss(Gaussian)')
end
%% Ensure stats match up
%V=-(k)^2*((RI+n_imm).^2-n_imm^2);
sdev=std(RI,[],'all');
avg=mean(RI,'all');
end
end
