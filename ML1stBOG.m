function evol=ML1stBOG(fxx,fyy,lambda,n_imm,dz,V,U_in,opt)
% This function computes the propagation of an initial field U_in through a
% scattering potential V using a first Born approximation for each slice of V. fxx,fyy are the spatial frequency grid which have
% been ifftshifted. n_imm is the index of refraction of the background
% medium prop_phs is the propagation phase for the angular spectrum. dz is
% the size of the pixel in the direction of propagation in microns.
prop_crop           = (fxx.^2 + fyy.^2 > (n_imm/lambda)^2==0);
prop_phs= 1i*2*pi*sqrt((n_imm/lambda)^2-(fxx.^2+fyy.^2));
prop_phs=prop_phs.*prop_crop;
prop=@(z) exp(prop_phs*z);
Mask=(n_imm/lambda)^2>1.01*((fxx.^2+fyy.^2));
U=U_in; % Initial Field
switch opt
    case 'out'
%S=zeros(size(V));
for i=1:size(V,3)
    S=U;
    U=ifft2(prop(dz).*(fft2(U))); % Prop dz to next layer
    %Us=ifft2(AngGreens(dz).*fftshift(fft2(U.*obj(:,:,i).*dz))); % Scattered field of nth layer
    Us=ifft2((fft2(S.*V(:,:,i))).*(AngGreensn(Mask,fxx,fyy,lambda,n_imm,prop_phs,dz)))*dz; % Scattered field of nth layer
    %Us=conv_fft2(Greens(dz),U.*obj(:,:,i).*dz,'same');
    U=U+Us;
    %S(:,:,i)=U;
end
evol=S;

    case 'Vol'
        S=zeros(size(V));
for i=1:size(V,3)
    S(:,:,i)=U;
    U=ifft2(prop(dz).*(fft2(U))); % Prop dz to next layer
    %Us=ifft2(AngGreens(dz).*fftshift(fft2(U.*obj(:,:,i).*dz))); % Scattered field of nth layer
    Us=ifft2((fft2(S(:,:,i).*V(:,:,i))).*(AngGreensn(Mask,fxx,fyy,lambda,n_imm,prop_phs,dz)))*dz; % Scattered field of nth layer
    %Us=conv_fft2(Greens(dz),U.*obj(:,:,i).*dz,'same');
    U=U+Us;
    %S(:,:,i)=U;
end
evol=S;


end




function AG=AngGreensn(Mask,fxx,fyy,lambda,n_imm,prop_phs,z)

AG= ((-1i.*exp(prop_phs.*z)./(4.*pi.*sqrt((n_imm/lambda)^2-(fxx.^2+fyy.^2))))).*Mask;
AG(isnan(AG)==1)=0;
end
end