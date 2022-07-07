function [theta1, theta2]=thetafctrhodiff_2011(rho,varargin)
% calcul selon Atkinson 2005 en supposant n reel 
% Retourne la valeur de theta (2 solutions) pour un rho donn�.
% 
% 
% EXEMPLE :
% 
% [theta1, theta2]=thetafctrho(rho,indice)
%       avec :  rho, degr� de polarisation pour lequel on recherche les theta
%               indice, indice du mat�riau (peut �tre r��l ou complexe)
%                       indice = n;    ou      indice = n*(1+ik);
%                                   n : indice du mat�riau
%                                   k : coeff. d'extinction du mat�riau
%               theta1, angle (en rad), 1� solution correspondant � rho
%               theta2, angle (en rad), 2� solution correspondant � rho
% 
% [theta1, theta2]=thetafctrho(rho,n,k)
%       avec :  rho, degr� de polarisation pour lequel on recherche les theta
%               n : indice du mat�riau
%               k : coeff. d'extinction du mat�riau
%               theta1, angle (en rad), 1� solution correspondant � rho
%               theta2, angle (en rad), 2� solution correspondant � rho
% 
% 
% 
% 
% 
% +++++++++++++++++
% + INDICE REEL : +
% +++++++++++++++++
% 
%   PAS ENCORE IMPLEMENTE
%   L'�quation pour un indice r�el (n) est la suivante Atkinson (2005) :
%  RhoDif= ((n-1/n).^2.*sin(theta).^2) ./ (2 + 2*n^2 - (n+1/n)^2*sin(theta).^2  + 4.*cos(theta).*sqrt(n^2-sin(theta).^2)  ) ;
% 
%   Soit le syst�me � r�soudre :
%   ----------------------------
%       syms a b c d e x
%  x=sin(theta) 
%   A=r*(n+1/n).^2+(n-1/n).^2;
%   B=-(2*r+2*r*n.^2);
%    a=(16*r*r-A.^2);
%   b=16*r*(1+n.^2)+2*A*B;
%   c=16*r^2*n.^2-B^2;
% 
%   S=a*x^4-b*x^2+c
%   solve(S) 

%       Les 4 solutions du syst�me pr�c�dent sont :
% 
% s1= -((2*n^2*r - 10*n^4*r - 10*n^6*r + 2*n^8*r + 2*n^2*r^2 + 6*n^4*r^2 + 6*n^6*r^2 + 2*n^8*r^2 - 4*n^3*r*(n^8*r^2 - n^8 + 8*n^6*r^2 - 8*n^6*r - 2*n^4*r^2 - 16*n^4*r + 18*n^4 + 8*n^2*r^2 - 8*n^2*r + r^2 - 1)^(1/2))/(n^8*r^2 + 2*n^8*r + n^8 + 4*n^6*r^2 - 4*n^6 - 10*n^4*r^2 - 4*n^4*r + 6*n^4 + 4*n^2*r^2 - 4*n^2 + r^2 + 2*r + 1))^(1/2)
% s2= -((2*n^2*r - 10*n^4*r - 10*n^6*r + 2*n^8*r + 2*n^2*r^2 + 6*n^4*r^2 + 6*n^6*r^2 + 2*n^8*r^2 + 4*n^3*r*(n^8*r^2 - n^8 + 8*n^6*r^2 - 8*n^6*r - 2*n^4*r^2 - 16*n^4*r + 18*n^4 + 8*n^2*r^2 - 8*n^2*r + r^2 - 1)^(1/2))/(n^8*r^2 + 2*n^8*r + n^8 + 4*n^6*r^2 - 4*n^6 - 10*n^4*r^2 - 4*n^4*r + 6*n^4 + 4*n^2*r^2 - 4*n^2 + r^2 + 2*r + 1))^(1/2)
% s3=  ((2*n^2*r - 10*n^4*r - 10*n^6*r + 2*n^8*r + 2*n^2*r^2 + 6*n^4*r^2 + 6*n^6*r^2 + 2*n^8*r^2 - 4*n^3*r*(n^8*r^2 - n^8 + 8*n^6*r^2 - 8*n^6*r - 2*n^4*r^2 - 16*n^4*r + 18*n^4 + 8*n^2*r^2 - 8*n^2*r + r^2 - 1)^(1/2))/(n^8*r^2 + 2*n^8*r + n^8 + 4*n^6*r^2 - 4*n^6 - 10*n^4*r^2 - 4*n^4*r + 6*n^4 + 4*n^2*r^2 - 4*n^2 + r^2 + 2*r + 1))^(1/2)
% s4=  ((2*n^2*r - 10*n^4*r - 10*n^6*r + 2*n^8*r + 2*n^2*r^2 + 6*n^4*r^2 + 6*n^6*r^2 + 2*n^8*r^2 + 4*n^3*r*(n^8*r^2 - n^8 + 8*n^6*r^2 - 8*n^6*r - 2*n^4*r^2 - 16*n^4*r + 18*n^4 + 8*n^2*r^2 - 8*n^2*r + r^2 - 1)^(1/2))/(n^8*r^2 + 2*n^8*r + n^8 + 4*n^6*r^2 - 4*n^6 - 10*n^4*r^2 - 4*n^4*r + 6*n^4 + 4*n^2*r^2 - 4*n^2 + r^2 + 2*r + 1))^(1/2)
%  
%       D'o� :
%           theta1 = asin(  s1  ); 
%           theta2 = asin(  s2  ); 
%           theta3 = asin(  s3  ); 
%           theta4 = asin(  s4  ); 
%        
%               NB : on garde ??
% 
% +++++++++++++++++++++
% + INDICE COMPLEXE : +
% +++++++++++++++++++++
% 
%   L'�quation pour un indice complexe (nc=n*(-1+ik)) est la suivante :
%   rho  = (2.*real(nc).*tan(theta).*sin(theta))./ (tan(theta).^2.*sin(theta).^2 + abs(nc).^2);
%   ou
%   rho  = (2.*n.*tan(theta).*sin(theta))       ./ (tan(theta).^2.*sin(theta).^2 + n.^2*(1+k.^2));
% 
%   Soit le systeme � r�soudre : 
%   ----------------------------
% 
%   A.^2.*rho-2.*n.*A+n.^2.*rho.*(1+k.^2) = 0
%           avec A = tan(theta).*sin(theta)
% 
%       Solutions de ce syst�me :
%           A1 = (2.*n - sqrt(4.*n^2.*(1-rho.^2.*(1+k^2)))) ./ (2.*rho);
%           A2 = (2.*n + sqrt(4.*n^2.*(1-rho.^2.*(1+k^2)))) ./ (2.*rho);
% 
%           A1 = tan(theta).*sin(theta)
%           ==> cos(theta).^2+A1.*cos(theta)-1=0
%               Solution de ce syst�me :
%                   theta1 = acos( (-A1+sqrt(A1.^2+4)) ./ 2 );
%                   theta2 = acos( (-A1-sqrt(A1.^2+4)) ./ 2 );
% 
%           et
% 
%           A2 = tan(theta).*sin(theta)
%           ==> cos(theta).^2+A2.*cos(theta)-1=0
%               Solution de ce syst�me :
%                   theta3 = acos( (-A2+sqrt(A2.^2+4)) ./ 2 );
%                   theta4 = acos( (-A2-sqrt(A2.^2+4)) ./ 2 );
% 
% On obtient 4 solutions, dont 2 r�elles et 2 complexes
% 


%if nargin==2 % indice reel 
    indice = cell2mat(varargin(1));
    if ( isreal(indice) )
        % si l'indice n est r��l
        n = indice;
        %disp('L''indice est r�el'),
        % cas o� rho ~= 0
    else
        n=real(indice); 
    end
        index = find(rho~=0);
        r=rho(index);
        % Les 4 solutions du syst�me pr�c�dent sont :
        
s1b = -((2^(1/2)*n*((r.*(r - 2*n*(n^8*r.^2 - n^8 + 8*n^6*r.^2 - 8*n^6*r - 2*n^4*r.^2 - 16*n^4*r + 18*n^4 + 8*n^2*r.^2 - 8*n^2*r + r.^2 - 1).^(1/2) + 3*n^2*r + 3*n^4*r + n^6*r - 5*n^2 - 5*n^4 + n^6 + 1))./((r + 1).*(r + 6*n^2*r + n^4*r - 2*n^2 + n^4 + 1))) .^(1/2))./(n^2 - 1));
 s1=sqrt(abs(s1b));
s2b= -((2^(1/2)*n*((r.*(r + 2*n*(n^8*r.^2 - n^8 + 8*n^6*r.^2 - 8*n^6*r - 2*n^4*r.^2 - 16*n^4*r + 18*n^4 + 8*n^2*r.^2 - 8*n^2*r + r.^2 - 1).^(1/2) + 3*n^2*r + 3*n^4*r + n^6*r - 5*n^2 - 5*n^4 + n^6 + 1))./((r + 1).*(r + 6*n^2*r + n^4*r - 2*n^2 + n^4 + 1)))  .^(1/2))./(n^2 - 1));
 s2=sqrt(abs(s2b));
s3b = ((2^(1/2)*n*((r.*(r - 2*n*(n^8*r.^2 - n^8 + 8*n^6*r.^2 - 8*n^6*r - 2*n^4*r.^2 - 16*n^4*r + 18*n^4 + 8*n^2*r.^2 - 8*n^2*r + r.^2 - 1).^(1/2) + 3*n^2*r + 3*n^4*r + n^6*r - 5*n^2 - 5*n^4 + n^6 + 1))./((r + 1).*(r + 6*n^2*r + n^4*r - 2*n^2 + n^4 + 1))).^(1/2))./(n^2 - 1));
s3=sqrt(abs(s3b)); 
% s3=asin(s3b);
s4b = ((2^(1/2)*n*((r.*(r + 2*n*(n^8*r.^2 - n^8 + 8*n^6*r.^2 - 8*n^6*r - 2*n^4*r.^2 - 16*n^4*r + 18*n^4 + 8*n^2*r.^2 - 8*n^2*r + r.^2 - 1).^(1/2) + 3*n^2*r + 3*n^4*r + n^6*r - 5*n^2 - 5*n^4 + n^6 + 1))./((r + 1).*(r + 6*n^2*r + n^4*r - 2*n^2 + n^4 + 1))).^(1/2))./(n^2 - 1));
s4=sqrt(abs(s4b)); 
% s4=asin(s4b);

 
% s1b= ((2*n^2*r - 10*n^4*r - 10*n^6*r + 2*n^8*r + 2*n^2*r.^2 + 6*n^4*r.^2 + 6*n^6*r.^2 + 2*n^8*r.^2 - 4*n^3*r.*(n^8*r.^2 - n^8 + 8*n^6*r.^2 - 8*n^6*r - 2*n^4*r.^2 - 16*n^4*r + 18*n^4 + 8*n^2*r.^2 - 8*n^2*r + r.^2 - 1).^(1/2))./(n^8*r.^2 + 2*n^8*r + n^8 + 4*n^6*r.^2 - 4*n^6 - 10*n^4*r.^2 - 4*n^4*r + 6*n^4 + 4*n^2*r.^2 - 4*n^2 + r.^2 + 2*r + 1));
% % s1=sqrt(abs(s1b));
% s1=-asin(s1b);
% s2b= ((2*n^2*r - 10*n^4*r - 10*n^6*r + 2*n^8*r + 2*n^2*r.^2 + 6*n^4*r.^2 + 6*n^6*r.^2 + 2*n^8*r.^2 + 4*n^3*r.*(n^8*r.^2 - n^8 + 8*n^6*r.^2 - 8*n^6*r - 2*n^4*r.^2 - 16*n^4*r + 18*n^4 + 8*n^2*r.^2 - 8*n^2*r + r.^2 - 1).^(1/2))./(n^8*r.^2 + 2*n^8*r + n^8 + 4*n^6*r.^2 - 4*n^6 - 10*n^4*r.^2 - 4*n^4*r + 6*n^4 + 4*n^2*r.^2 - 4*n^2 + r.^2 + 2*r + 1));
% % s2=sqrt(abs(s2b));
%  s2=-asin(s2b);
% s3b=  ((2*n^2*r - 10*n^4*r - 10*n^6*r + 2*n^8*r + 2*n^2*r.^2 + 6*n^4*r.^2 + 6*n^6*r.^2 + 2*n^8*r.^2 - 4*n^3*r.*(n^8*r.^2 - n^8 + 8*n^6*r.^2 - 8*n^6*r - 2*n^4*r.^2 - 16*n^4*r + 18*n^4 + 8*n^2*r.^2 - 8*n^2*r + r.^2 - 1).^(1/2))./(n^8*r.^2 + 2*n^8*r + n^8 + 4*n^6*r.^2 - 4*n^6 - 10*n^4*r.^2 - 4*n^4*r + 6*n^4 + 4*n^2*r.^2 - 4*n^2 + r.^2 + 2*r + 1));
% s3=sqrt(abs(s3b));
% s4b=  ((2*n^2*r - 10*n^4*r - 10*n^6*r + 2*n^8*r + 2*n^2*r.^2 + 6*n^4*r.^2 + 6*n^6*r.^2 + 2*n^8*r.^2 + 4*n^3*r.*(n^8*r.^2 - n^8 + 8*n^6*r.^2 - 8*n^6*r - 2*n^4*r.^2 - 16*n^4*r + 18*n^4 + 8*n^2*r.^2 - 8*n^2*r + r.^2 - 1).^(1/2))./(n^8*r.^2 + 2*n^8*r + n^8 + 4*n^6*r.^2 - 4*n^6 - 10*n^4*r.^2 - 4*n^4*r + 6*n^4 + 4*n^2*r.^2 - 4*n^2 + r.^2 + 2*r + 1));
% s4=sqrt(abs(s4b));

        % D'o� :
        theta1(index) = asin( s3);   % r��l au poui�me pr�s
        theta2(index) = asin( s4 );   % r��l au poui�me pr�s
        theta1(index) = real(theta1(index));   % retrait du poui�me imaginaire
        theta2(index) = real(theta2(index));   % retrait du poui�me imaginaire
        % cas o� rho == 0 ==> theta1=0 et theta2=0
        index = find(rho==0);
        theta1(index)=0;
        theta2(index)=0;
        
%     else
%         % si l'indice n est complexe;  avoir dans ce cas ! 
%         nc = indice; % nc = n.*(1+i.*k);
%         n = real(nc);
%         k = imag(nc)/n;
%         %disp('L''indice est complexe'),
%         % cas o� rho ~= 0
%         index = find(rho~=0);
%         A1(index) = (2.*n - sqrt(4.*n^2.*(1-rho(index).^2.*(1+k^2)))) ./ (2.*rho(index));
%         A2(index) = (2.*n + sqrt(4.*n^2.*(1-rho(index).^2.*(1+k^2)))) ./ (2.*rho(index));
%         % cas o� rho == 0 ==> A1=0 et A2=0
%         index = find(rho==0);
%         A1(index)=0;
%         A2(index)=0;
% %         theta1 = acos( (-A1+sqrt(A1.^2+4)) ./ 2 );
% %         theta2 = acos( (-A1-sqrt(A1.^2+4)) ./ 2 ); % complexe
% %         theta3 = acos( (-A2+sqrt(A2.^2+4)) ./ 2 );
% %         theta4 = acos( (-A2-sqrt(A2.^2+4)) ./ 2 ); % complexe
%         cossol1 = ((-A1+sqrt(A1.^2+4)) ./ 2); % -1 <= cossol1 <= 1
%         cossol2 = ((-A1-sqrt(A1.^2+4)) ./ 2); % en dehors des bornes de acos
%         cossol3 = ((-A2+sqrt(A2.^2+4)) ./ 2); % -1 <= cossol1 <= 1
%         cossol4 = ((-A2-sqrt(A2.^2+4)) ./ 2); % en dehors des bornes de acos
%         % test si  -1 <= cossolx <= 1
%         test1 = sum( cossol1>=-1 & cossol1<=1 ) == max(size(cossol1(~isnan(cossol1))));
%         test2 = sum( cossol2>=-1 & cossol2<=1 ) == max(size(cossol2(~isnan(cossol2))));
%         test3 = sum( cossol3>=-1 & cossol3<=1 ) == max(size(cossol3(~isnan(cossol3))));
%         test4 = sum( cossol4>=-1 & cossol4<=1 ) == max(size(cossol4(~isnan(cossol4))));
%         if test1
%             theta1 = acos( cossol1 );
%             solution1 = theta1;
%         end
%         if test2
%             theta2 = acos( cossol2 ); % complexe
%             solution1 = theta2;
%         end
%         if test3
%             theta3 = acos( cossol3 );
%             solution2 = theta3;
%         end
%         if test4
%             theta4 = acos( cossol4 ); % complexe
%             solution2 = theta4;
%         end
%         theta1 = solution1;
%         theta2 = solution2;
%     end
% elseif nargin==3
%     n = cell2mat(varargin(1));
%     k = cell2mat(varargin(2));
%         nc = n.*(1+i.*k);
%         %disp('L''indice est complexe'),
%         % cas o� rho ~= 0
%         index = find(rho~=0);
%         A1(index) = (2.*n - sqrt(4.*n^2.*(1-rho(index).^2.*(1+k^2)))) ./ (2.*rho(index));
%         A2(index) = (2.*n + sqrt(4.*n^2.*(1-rho(index).^2.*(1+k^2)))) ./ (2.*rho(index));
%         % cas o� rho == 0 ==> A1=0 et A2=0
%         index = find(rho==0);
%         A1(index)=0;
%         A2(index)=0;
% %         theta1 = acos( (-A1+sqrt(A1.^2+4)) ./ 2 );
% %         theta2 = acos( (-A1-sqrt(A1.^2+4)) ./ 2 ); % complexe
% %         theta3 = acos( (-A2+sqrt(A2.^2+4)) ./ 2 );
% %         theta4 = acos( (-A2-sqrt(A2.^2+4)) ./ 2 ); % complexe
%         cossol1 = ((-A1+sqrt(A1.^2+4)) ./ 2); % -1 <= cossol1 <= 1
%         cossol2 = ((-A1-sqrt(A1.^2+4)) ./ 2); % en dehors des bornes de acos
%         cossol3 = ((-A2+sqrt(A2.^2+4)) ./ 2); % -1 <= cossol1 <= 1
%         cossol4 = ((-A2-sqrt(A2.^2+4)) ./ 2); % en dehors des bornes de acos
%         % test si  -1 <= cossolx <= 1
%         test1 = sum( cossol1>=-1 & cossol1<=1 ) == max(size(cossol1(~isnan(cossol1))));
%         test2 = sum( cossol2>=-1 & cossol2<=1 ) == max(size(cossol2(~isnan(cossol2))));
%         test3 = sum( cossol3>=-1 & cossol3<=1 ) == max(size(cossol3(~isnan(cossol3))));
%         test4 = sum( cossol4>=-1 & cossol4<=1 ) == max(size(cossol4(~isnan(cossol4))));
%         if test1
%             theta1 = acos( cossol1 );
%             solution1 = theta1;
%         end
%         if test2
%             theta2 = acos( cossol2 ); % complexe
%             solution1 = theta2;
%         end
%         if test3
%             theta3 = acos( cossol3 );
%             solution2 = theta3;
%         end
%         if test4
%             theta4 = acos( cossol4 ); % complexe
%             solution2 = theta4;
%         end
%         theta1 = solution1;
%         theta2 = solution2;
% end

[m,n] = size(rho);
theta1 = reshape(theta1,m,n);
theta2 = reshape(theta2,m,n);

% suppression du pouieme d'imaginaire
theta1 = real(theta1);
theta2 = real(theta2);

