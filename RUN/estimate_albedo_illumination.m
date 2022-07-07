function [albedo,I_d,slant,tilt] = estimate_albedo_illumination (dif)
% E: the image
% normalizing the image to have maximum of one
%E = E ./ max(E(:));
% compute the average of the image brightness
Mu1 = mean(dif(:));
% compute the average of the image brightness square
Mu2 = mean(mean(dif.^2));
 
% now lets compute the image's spatial gradient in x and y directions
[Ex,Ey] = gradient(dif);
% %imshow(Ex)
% %imshow(Ey)
% figure(1),
% subplot(1,2,1), imshow(Ex)
% subplot(1,2,2), imshow(Ey)
 
% normalize the gradients to be unit vectors
Exy = sqrt(Ex.^2 + Ey.^2);
nEx = Ex ./(Exy + eps); % to avoid dividing by zero
nEy = Ey ./(Exy + eps);
% % imshow(nEx)
% % imshow(nEy)
% figure(2),
% subplot(1,2,1), imshow(nEx)
% subplot(1,2,2), imshow(nEy)
 
% computing the average of the normalized gradients
avgEx = mean(nEx(:));
avgEy = mean(nEy(:));
 
% now lets estimate the surface albedo
gamma = sqrt((6 *(pi^2)* Mu2) - (48 * (Mu1^2)));
albedo = gamma/pi;
% estimating the slant
slant = acos((4*Mu1)/gamma);
% estimating the tilt
tilt = atan(avgEy/avgEx);
if tilt < 0
    tilt = tilt + pi;
end
% the illumination direction will be:
I_d = [cos(tilt)*sin(slant) sin(tilt)*sin(slant) cos(slant)];
end
