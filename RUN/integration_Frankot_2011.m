function z = integration_Frankot(Ip,Iq,dx,dy)

z = integration_Frankot_mathias(Ip,Iq,dx,dy);
z = morel(Ip,Iq,dx,dy);
needlemap(:,:,1)=-Ip;
needlemap(:,:,2)=-Iq;
needlemap(:,:,3)=ones(size(Ip));
z = FrankotChellappa(needlemap);
% display needle map


function Iz = integration_Frankot_mathias(Ip,Iq,dx,dy)
%dx=-dx;

[lig,col] = size(Ip);

%% prepa des axes
for i1=1:lig
    for j1=1:col
        u(i1,j1) = (j1-1)-col./2;
        v(i1,j1) = (i1-1)-lig./2;
    end
end
u = u./col.*(2*pi)./dx;
v = v./lig.*(2*pi)./dy;
u = fftshift(u);
v = fftshift(v);
        %figure, hold on, imshow(u,[]), colorbar, title('u'),
        %figure, hold on, imshow(v,[]), colorbar, title('v'),
%% transfo de fourier du gradient
tfp = fft2(Ip);
tfq = fft2(Iq);
%         figure, hold on, imshow(real(tfp),[]), colorbar, title('re(tfp)'),
%         figure, hold on, imshow(imag(tfp),[]), colorbar, title('im(tfp)'),
%         figure, hold on, imshow(log(abs(tfp)+1),[]), colorbar, title('tfp'),
%         figure, hold on, imshow(log(abs(tfq)+1),[]), colorbar, title('tfq'),

%% integration
A  = -i.*( tfp.*u + tfq.*v );
B  = (u.^2 + v.^2);
%disp('nb de 0 de B'), max(size(find(B==0)))
B(find(B==0)) = eps;%1;%0.0000001; % ne pas diviser par 0
tmp = A ./ B;
%         figure, hold on, imshow(log(abs(tmp)+1),[]), colorbar, title('tmp'),
%%
tmp2 = real(ifft2(tmp));
%         figure, hold on, imshow(tmp2,[]), colorbar, title('tmp2'),
%         figure, contour(tmp2,15), colormap(copper), title('tmp2-contour'),

%% pentes
[lig,col] = size(tfq);
tmpp = dx*real(tfp(1,1))/(lig*col);
tmpq = dy*real(tfq(1,1))/(lig*col);
for i1=1:lig
    for j1=1:col
        pentes(i1,j1) = (j1-1)*tmpp + (i1-1)*tmpq;
    end
end
%         figure, hold on, imshow(pentes,[]), colorbar, title('pentes'),
%         figure, contour(pentes,15), colormap(copper), title('pentes-contour'),
%% z
Iz = pentes + tmp2;
%         figure, hold on, imshow(Iz,[]), colorbar, title('Iz'),
%         figure, contour(Iz,15), colormap(copper), title('Iz-contour'),
        
        
        
        
 

function z=morel(p,q,varargin)       
% MOREL  - Generates integrable surface from gradients
%
% method of Frankot-Chellappa
% using FFT
%
% Usage:      z = morel(p,q,dx,dy)
%
% Arguments:  p,  - 2D matrices specifying a grid of gradients of z
%             q     with respect to x and y.
%
% Options :   dx, - step between pixels
%               dy
%
% Returns:    z      - Inferred surface heights.

% Reference:
%
% O. Morel

% Olivier Morel
% Le2i
%
% January 2005

[nbl,nbc]=size(p);
dx=1;
dy=1;

if length(varargin)>= 1
   dx=varargin{1};
   dy=dx;
end
if length (varargin) == 2
   dy=varargin{2};
end

% Fourier transform of the gradient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tf_p=fft2(p);
tf_q=fft2(q);

%Pr�paration des axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=fftshift(((0:nbc-1)-nbc/2)/nbc)/dx*2*pi;
v=fftshift(((0:nbl-1)-nbl/2)/nbl)/dy*2*pi;
for k=1:nbl
   U(k,:)=u;
end
for k=1:nbc
   V(:,k)=v';
end

%int�gration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfZ=-i*(U.*tf_p+V.*tf_q)./(U.^2+V.^2+eps);

%constante pente
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pente_x=real(tf_p(1,1))/(nbc*nbl);
pente_y=real(tf_q(1,1))/(nbc*nbl);

%calcul des axes
x=0:dx:dx*(nbc-1);
y=0:dy:dy*(nbl-1);

for l=1:nbl
   for c=1:nbc
      const(l,c)=pente_x*x(c)+pente_y*y(l);
   end
end
z=real(ifft2(tfZ))+const;
