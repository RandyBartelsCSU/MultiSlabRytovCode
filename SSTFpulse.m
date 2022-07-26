function E_sstf=SSTFpulse(E0,E_w,Omega,omega0,NA,x,y,z,n_imm,alpha,F2)
 % This function outputs the field at a location z in space from the focus for a simultaneous
 % space time focussing gaussian beam for a certain angular frequency
 % Omega where Omega=omega-omega0
 
% E_w is the spectrum of the pulse (in w-w0) with initial spectral phase if
% desired

% E0 is the amplitude of the field

% alpha is the spatial chirp parameter

% F2 is the focal length of 
c=(3e8)*1e6*1e-15;                                 % speed of light microns per fs
k=(Omega/c)*n_imm;
k0=omega0/c;
lambda=2*pi*c/Omega;

% Gaussian Beam Parameters
        w0=.61*lambda/asin(NA/n_imm);              % Beam waist function of wavelength                                          
        z0=pi*n_imm*w0^2/lambda;                   % Rayleigh range
        R=@(z0,z) z.*(1+(z0./z).^2);               % Radius of curvature of wavefront
        phi=@(z,z0) atan(z/z0);                    % Axial phase
        w=@(z0,w0,z) sqrt(w0^2*(1+(z/z0).^2));     % Beam Radius
        
        
        
        %@(x,y,z,omega,omega0,z0,w0,alpha,F2)
Phi_sstf=k.*z-(k/2)*((alpha.*Omega./F2).^2).*z-k.*(alpha.*Omega./F2).*x+k.*(((x+(alpha.*Omega./F2).*z).^2)./(2.*R(z0,z)))+(k.*y.^2)./(2.*R(z0,z))-phi(z,z0);
E_sstf=E0*(w0./w(z0,w0,z)).*E_w.*exp(-(((x+(alpha.*Omega./F2).*z).^2+y.^2)./(w(z0,w0,z).^2))).*exp(1i.*Phi_sstf);
end
