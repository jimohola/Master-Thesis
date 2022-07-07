clc;
close all;
clear all;


I0 = imread('I0.jpg');
I0= imcrop(I0,[224 393 1139 967]);
I0 = imgaussfilt(I0,4);

I45 = imread('I45.jpg');
I45= imcrop(I45,[224 393 1139 967]);
I45 = imgaussfilt(I45,4);

I90 = imread('I90.jpg');
I90= imcrop(I90,[224 393 1139 967]);
I90 = imgaussfilt(I90,4);

I135  = imread('I135.jpg');
I135= imcrop(I135,[224 393 1139 967]);
I135 = imgaussfilt(I135,4);

% [I0 rec] = imcrop(imread('I0.jpg'));
% I45 =  imcrop(imread('I45.jpg'), rec);
% I90 =  imcrop(imread('I90.jpg'), rec);
% I135 =  imcrop(imread('I135.jpg'), rec);

I0 = double(rgb2gray(I0));
I45 = double(rgb2gray(I45));
I90 = double(rgb2gray(I90));
I135  = double(rgb2gray(I135));

figure;
subplot(221), imagesc(I0); title('I0');
subplot(222), imagesc(I45); title('I45');
subplot(223), imagesc(I90); title('I90');
subplot(224), imagesc(I135);  title('I135');


%% 

T = (I0 + I90);
phi = (1/2)*atan2((T - 2*I45), (I90-I0)) + pi/2;
rho = abs( (I90 - I0)./((I90 + I0).*cos(2*phi)) );

figure;
subplot(1,3,1);imagesc(T); stitle=sprintf('T'); title(stitle);
subplot(1,3,2);imagesc((rho)); stitle=sprintf('rho'); title(stitle);
subplot(1,3,3);imagesc(phi); stitle=sprintf('phi'); title(stitle);

%% Gpolar
n=1.4;
[theta1, theta2]=thetafctrhodiff_2011(rho,n);

PHI= phi;
ind=  I90 < I0 ;
PHI(ind)=phi(ind)+pi/2;

ind =  I45 < I0 ;
PHI(ind)=phi(ind)+pi;

figure;
Gx = tan(theta2).*sin(PHI);  % x component for normal vector
Gy = tan(theta2).*cos(PHI); % y component for normal vector
imshowpair(Gx,Gy,'montage'); stitle=sprintf('Gpolar'); title(stitle);


%% ICA algorithm (Mutual Information)

[row,col] = size(I0);
X = [I0(:)' ; I45(:)';  I90(:)'  ; I135(:)'];
% X = [I135(:)' ; I90(:)';  I45(:)'  ; I0(:)' ];

% from svd of X, find W such that X=UcWW-1 Dc Vc'
[U D V]= svds(X,2);

%first column of W
ab=U'*ones(length(U),1); % coloumn contain a and b from matrix W
% ab=U'*[1; 1/2; 0; 1/2] ; % coloumn contain a and b from matrix W

a =ab(1,1); b = ab(2,1); %ab value of a and b

alpha =atan (b/a); %find alpha from value a and b
% alpha=atan2(b,a)

R1 = abs ( a/cos (alpha) ); %positive value of r1 

k=1; %represent array for I, mutual information

%%
mutu_info =0; I=[]; 
Beta=0:0.01:pi;
for beta = 0:0.01:pi
    R2 =abs( 1/(R1*sin(beta-alpha))) ;
    a = R1*cos(alpha);
    b = R1*sin(alpha);
    c = R2*cos(beta);
    d = R2*sin(beta);
    
    W = [a c ; 
             b d];
         
    S = inv(W) * D *V';
    diff = abs ( S (1,:));
    spclr = abs ( S(2,:));

    mutu_info = joint_h2011 (diff,spclr,row,col);

    I(k) = mutu_info;
    k = k+1;
end

figure;
plot(0:0.01:pi, I), 
colorbar, title('Mutual information');

[m ind] = min(I);
beta2=Beta(ind);

%%
R2 =abs( 1/(R1*sin(beta2-alpha))) ;
    a = R1*cos(alpha);
    b = R1*sin(alpha);
    c = R2*cos(beta2);
    d = R2*sin(beta2);
    
    W2 = [a c ; b d];

    S = inv(W2) * D *V';
    
    diff = abs ( S (1,:));
    spclr = abs ( S(2,:));
    
    diff2  = zeros (row,col);
    spclr2 = zeros (row,col);

    diff2  = reshape (diff,row,col);
    spclr2 = reshape (spclr,row,col);

    dif = abs( diff2); 
    spec =  abs( spclr2);

figure; 
imshow(dif,[]), colorbar, title('(a)  Diffuse');
figure;
imshow(spec,[]), colorbar, title('(b)  specular');


% imshowpair(dif,spec,'montage'); stitle=sprintf('Separated Mixed Reflection'); title(stitle);

%% Gradient field 

% first compute the surface albedo and illumination direction
[albedo,I_d,slant,tilt] = estimate_albedo_illumination (dif);

% initializations
[M,N] = size(dif);

% surface normals
p = zeros(M,N);
q = zeros(M,N);


% the surface
Z = zeros(M,N);

% surface derivatives in x and y directions
Z_x = zeros(M,N);
Z_y = zeros(M,N);

% maximum number of iterations
maxIter = 200;

% the normalized illumination direction
ix = cos(tilt) * tan(slant);
iy = sin(tilt) * tan(slant);
 
for k = 1 : maxIter
    % using the illumination direction and the currently estimate surface normals, compute the corresponding reflectance map.
    R =(cos(slant) + p .* cos(tilt)*sin(slant)+ q .* ...
        sin(tilt)*sin(slant))./sqrt(1 + p.^2 + q.^2);
    
    % at each iteration, make sure that the reflectance map is positive at each pixel, set negative values to zero.
    R = max(0,R);
    
    % compute our function f which is the deviation of the computed reflectance map from the original image ...
    f = dif - R;
    
    % compute the derivative of f with respect to our surface Z
    df_dZ =(p+q).*(ix*p + iy*q + 1)./(sqrt((1 + p.^2 + q.^2).^3)* ...
        sqrt(1 + ix^2 + iy^2))-(ix+iy)./(sqrt(1 + p.^2 + q.^2)* ...
        sqrt(1 + ix^2 + iy^2));
    
    % update our surface
    Z = Z - f./(df_dZ + eps); % to avoid dividing by zero
    
    % compute the surface derivatives with respect to x and y
    Z_x(2:M,:) = Z(1:M-1,:);
    Z_y(:,2:N) = Z(:,1:N-1);
    
    % using the updated surface, compute new surface normals
    p = Z - Z_x;
    q = Z - Z_y;
end


% smoothing the recovered surface
Z = medfilt2((abs(Z)),[21 21]);
[Zx Zy] = imgradient(Z);
figure;
imshowpair(Zx,Zy,'montage'); stitle=sprintf('Gdepth'); title(stitle);
% imshowpair(p,q,'montage');

% visualizing the result
figure;
surfl(Z);
shading interp;
colormap gray(256);
lighting phong;

% view gradients  
% figure
% subplot(211)
% imagesc( ((-1.5e20 < p)&(p < 1.5e20)).*p );
% subplot(212)
% imagesc( ((-1.5e20 < q)&(q < 1.5e20)).*q );

% subsampling factor to view normals 
subs=8;
p=p(1:subs:M,1:subs:N);
Q=q(1:subs:M,1:subs:N);

R1=ones(round(M/subs),round(N/subs));


%% Correction

V1= vecnorm([ Zx(:)'-Gx(:)'; Zy(:)'-Gy(:)']);
V2= vecnorm([ Zx(:)'+Gx(:)'; Zy(:)'+Gy(:)']);

T1=V2 < V1 ; 
T1=reshape(T1, M,N);

CorPHI=PHI+T1*pi; 

% Gcorrect
figure;
Gcx = tan(theta2).*sin(CorPHI);  % x component for normal vector
Gcy= tan(theta2).*cos(CorPHI); % y component for normal vector
imshowpair(Gcx,Gcy,'montage'); stitle=sprintf('Gcorrect'); title(stitle);



%% calculate the difference between Gpolar and Gcor

%extract the size of the gradients image
[X, Y]=size(Gcx);   % sine Gcx

Nx=[Gx(:)'-Zx(:)';  Gy(:)'-Zy(:)'; ones( X*Y )]; % a 3 by X*Y matrix
Norms=vecnorm(Nx);
INorms=reshape(Norms,X,Y); %
imshow(INorms);


%% z before correction 

% dx =0.025;
% dy = 0.025;
% z = integration_Frankot_2011(Gx,Gy,dx,dy);

nb_iter= 1000;
filter = [1 1 1 1 1;1 1 1 1 1;1 1 0 1 1;1 1 1 1 1;1 1 1 1 1];
z = relaxation(Gx,Gy,nb_iter,filter);
 
figure;
surfl(z);
shading interp;
colormap gray(256);
lighting phong;
stitle=sprintf('Surface recovery before Correction'); title(stitle);

%% z after correction 

% dx =0.025;
% dy = 0.025;
% z = integration_Frankot_2011(Gcx,Gcy,dx,dy);
% z = medfilt2((abs(z)),[21 21]);

nb_iter= 1000;
filter = [1 1 1 1 1;1 1 1 1 1;1 1 0 1 1;1 1 1 1 1;1 1 1 1 1];
z = relaxation(Gcx,Gcy,nb_iter,filter);


figure;
surfl(z);
shading interp;
colormap gray(256);
lighting phong;
stitle=sprintf('Surface recovery after Correction'); title(stitle);

