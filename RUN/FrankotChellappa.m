function [Z,PQ,Fourier] = FrankotChellappa(needlemap, varargin)
%FrankotChellappa v2
%
% Usage
% * |FrankotChellapa(/needlemap/, /options/)| Integrates the needlemap
%   into a surface using Frankot and Chellappa's method.
%
% Notes
% * This algorithm is based on the following seminal paper:
%   (1) R.T.Frankot and Z.Chellappa, "A method for enforcing integrability    
%   in shape from shading algorithms" IEEE Transactions in Pattern
%   Recognition and Machine Intelligence (1988) 10. pp. 439-451
% * The original author of this code is Mario Castelan, Department of
%   Computer Science.  University of York, Heslington, York, YO10 5DD, UK
% * Contains code by C. Thompson and S. Eddins, from The Mathworks, Roberto 
%   Fraile, from cs.york.ac.uk, and ASG.
% * Option |Transform="Cosine"| uses the Discrete Cosine Transform instead
%   of the Fast Fourier Transform.  This is the default.
% * Option |Transform="Fourier"| uses the Fast Fourier Transform.
% * Option |GradientFilter=/number/| filters the gradients that will be
%   taken into account.  The default value is 4 (As the basis functions to
%   minimize the equation proposed in (1) are the Fourier transform, special
%   care should be taken for values where a high slant is presented, that is
%   boundaries, these ones are set to zero).
% * Options |Delta=/number/| and |Miu=/number/| set regularising terms referring
%   to curvature and area constraints as proposed in (2).  Their default value is zero,
%   which corresponds to the algorithm in the original paper (1).
%   (2) T.Wei and R.Klette. "Height from gradients with surface curvature and area
%   constraints" (Technical report) CITR-TR-109, Tamaki Campus, University of 
%   Auckland, October 2001. 
% * A more complete explanation about implementing this method can be found in
%   (3) R.Klette,K.Schluns,A.Koschan "Computer Vision: Three-Dimensional Data from Images"  
%   Springer-Verlag Singapore Pte. Ltd. ISBN: 9813083719



% hardcoded values as advised in the documentation (these should become options)
option_GradientFilter = -1; %4;
disp('FC algo; no gradient filtering')
option_Delta = 0;
option_Miu = 0;

% default values
%option_Transform = 'Cosine';
option_Transform = 'Fourier';

%DeveloperParseOptions(varargin);

% get the gradient from the needlemap
PQ = convPQ(needlemap);

switch (option_Transform)
 case 'Fourier'
  
  if option_GradientFilter == -1 %  special case when one does not want the gradient field to be filtered.
    option_GradientFilter = max(max(PQ));
  end;   
  
  Y = size(PQ,1);  % the image is supposed to be even, for the antysimetr
  X = size(PQ,2);
  
  for y = 1:Y     % filtering the gradients
    for x = 1:X
      if  abs(PQ(y,x,1)) <=  option_GradientFilter & abs(PQ(y,x,2)) <=  option_GradientFilter
        P1(y,x) = PQ(y,x,1);
        Q1(y,x) = PQ(y,x,2);
        P2(y,x) = 0;
        Q2(y,x) = 0;            
      else
        P1(y,x) = 0;
        P2(y,x) = 0;
        Q1(y,x) = 0;
        Q2(y,x) = 0;
      end;    
    end;
  end;   
  
  
  FouP = fft2(P1); FouQ = fft2(Q1); % transforming the gradient fields to the frecuency domain
  P1 = real(FouP); Q1 = real(FouQ);
  P2 = imag(FouP); Q2 = imag(FouQ);
  
  H1 = zeros(Y,X);
  H2 = zeros(Y,X);
  
  for y = 1:Y % Here the integration process using the Fourier basis functions starts
    v = y-1;
    
    for x = 1:X
      u = x-1;
      
      if ~(v == 0 & u == 0)
        
        if u <= ((X/2)-1)
          AxWx = (2*pi*u)/X;
        else
          AxWx = ((2*pi*u)/X)-(2*pi);
        end;
        if v <= ((Y/2)-1)
          AyWy = (2*pi*v)/Y;
        else
          AyWy = ((2*pi*v)/Y)-(2*pi);
        end;
        
        %AxWx = sin((2 * pi * u)/(tam));  %this can also be used as an approximation of the differentiation
        %AyWy =  sin((2 * pi * v)/(tam)); %operator for the Fourier transform but only when the range of (u,v)
        %                                  is from -M/2 to (M/2)-1, assuming even M (not this case where 0 <= [u,v] >= M-1) 
        
        H1(y,x) = ((AxWx*P2(y,x)) + (AyWy*Q2(y,x)))/(((1+option_Delta)*(AxWx^2+AyWy^2))+(option_Miu*(AxWx^2+AyWy^2)^2));
        H2(y,x) = ((-AxWx*P1(y,x)) + (-AyWy*Q1(y,x)))/(((1+option_Delta)*(AxWx^2+AyWy^2))+(option_Miu*(AxWx^2+AyWy^2)^2));
        
      end;
    end;
  end;    
  
  H1(1,1) = 0;
  H2(1,1) = 0;
  
  Fourier = complex(H1,H2);
  
  Z = ifft2(Fourier);
  Z = -real(Z);             % The real part of the inverse transform contains the height map
%   Z = Z + max(max(abs(Z))); % normalization
%   Z = Z/max(max(Z));
  
  %  This part of the code can be activated if one wants to recover also the new
  %  -nearest integrable- gradient field [Zx,Zy] corresponding to the new recovered surface
  
  for y = 1:Y 
    v = y-1;
    for x = 1:X
      u = x-1;
      if ~(v == 0 & u == 0)
        if u <= ((X/2)-1)
          AxWx = (2*pi*u)/X;
        else
          AxWx = ((2*pi*u)/X)-(2*pi);
        end;
        if v <= ((Y/2)-1)
          AyWy = (2*pi*v)/Y;
        else
          AyWy = ((2*pi*v)/Y)-(2*pi);
        end;
        Zx(y,x) =  -j * AxWx * Fourier(y,x);
        Zy(y,x) =  -j * AyWy * Fourier(y,x);
      end;
    end;
  end;    
  Zx(1,1) = 0;
  Zy(1,1) = 0;
  
  PQ(:,:,1) = -real(ifft2(Zx));
  PQ(:,:,2) = -real(ifft2(Zy));
  
  % setting the background as zero if needed

%for y = 1:Y
%    for x = 1:X
%        if PQ(y,x,1) == 0 & PQ(y,x,2) == 0
%           H(y,x) = 0; 
%       end;
%end;
%end;

 case 'Cosine'
  
  Y = size(PQ,1);
  X = size(PQ,2);
  
  % filtering the gradient fields
  
  if option_GradientFilter ~= -1
    for y = 1:Y
      for x = 1:X
        if abs(PQ(y,x,1)) >  option_GradientFilter & abs(PQ(y,x,2)) >  option_GradientFilter        
          PQ(y,x,:) = 0; 
        end;
      end;
    end;
  end;
  
  
  % using the Discrete Cosine Transform to represent Zx and the Discrete Sine Transform DST o represent Zy
  
  tQx=dct(dst(PQ(:,:,1)')'); 
  tQy=dst(dct(PQ(:,:,2)')');
  
  % here the integration starts
  
  for y = 1:Y 
    AyWy = pi * y / Y; % the differentiation operator is approximated as pi * (u,v)/(M,N)
    for x = 1:X    
      AxWx = pi * x / X;
      
      c(y,x) =(-AxWx*tQx(y,x) - AyWy*tQy(y,x))/(((1+option_Delta)*(AxWx^2+AyWy^2))+(option_Miu*(AxWx^2+AyWy^2)^2)); 
      
    end;
  end;
  
  DC = c;
  Z = idct2(c);  % Height map delivered by the Inverse Discrete Cosine Transform
  
  Z = Z + max(max(abs(Z))); % normalizing
  Z = Z/max(max(Z));
  
  % setting the background as zero if needed
  
  %for y = 1:Y
  %    for x = 1:X
  %        if PQ(y,x,1) == 0 & PQ(y,x,2) == 0
  %           H(y,x) = 0; 
  %       end;
  %    end;
  %end;
  
 otherwise
  Print('FrankotChellappa doesn''t know the transform you are trying to use, please check the spelling');
end 



function PQ = convPQ(N)

tam1 = size(N,1);
tam2 = size(N,2);

PQ = zeros(tam1,tam2,2);
for x=1:tam1
    for y=1:tam2
        if N(x,y,3)==0
          PQ(x,y,1) = 0;
          PQ(x,y,2) = 0;
        else
          PQ(x,y,1) = (N(x,y,1)/N(x,y,3));
          PQ(x,y,2) = (N(x,y,2)/N(x,y,3));
        end;    
    end;
end;    


function b = dst(a,n)
%DST  Discrete sine transform.
%
%   Y = DST(X) returns the discrete sine transform of X.
%   The vector Y is the same size as X and contains the
%   discrete sine transform coefficients.
%
%   Y = DST(X,N) pads or truncates the vector X to length N 
%   before transforming.
%
%   If X is a matrix, the DST operation is applied to each
%   column.  This transform can be inverted using IDCT.
%
%   See also: FFT, IFFT, DCT, and IDCT.

%   Author(s): C. Thompson, 2-12-93
%              S. Eddins, 10-26-94, revised
%   Copyright (c) 1988-97 by The MathWorks, Inc.
%   $Revision: 1.14 $  $Date: 1997/02/06 21:52:37 $

%   Modified by ASG from DCT to DST.  1999/06/02

%   References: 
%   1) A. K. Jain, "Fundamentals of Digital Image
%      Processing", pp. 150-153.
%   2) Wallace, "The JPEG Still Picture Compression Standard",
%      Communications of the ACM, April 1991.

error(nargchk(1,2,nargin));

if min(size(a))==1 
    if size(a,2)>1
        do_trans = 1;
    else
        do_trans = 0;
    end
    a = a(:); 
else  
    do_trans = 0; 
end
if nargin==1,
  n = size(a,1);
end
m = size(a,2);

% Pad or truncate a if necessary
if size(a,1)<n,
  aa = zeros(n,m);
  aa(1:size(a,1),:) = a;
else
  aa = a(1:n,:);
end

if rem(n,2)==1 | ~isreal(a), % odd case
  % Form intermediate even-symmetric matrix.
  y = zeros(2*n,m);
  y(1:n,:) = aa;
  y(n+1:n+n,:) = -flipud(aa);

  % Perform FFT
  yy = fft(y);

  % Compute DST coefficients
  ww = (exp(-i*(0:n-1)*pi/(2*n))/sqrt(2*n)).';
  ww(1) = ww(1) / sqrt(2);
  b = ww(:,ones(1,m)).*yy(1:n,:);

else % even case

  % Re-order the elements of the columns of x
  y = [ aa(1:2:n,:); -aa(n:-2:2,:) ];

  % Compute weights to multiply DFT coefficients
  ww = 2*exp(-i*(0:n-1)'*pi/(2*n))/sqrt(2*n);
  ww(1) = ww(1) / sqrt(2);
  W = ww(:,ones(1,m));

  % Compute DST using equation (5.92) in Jain
  b = W .* fft(y);
end

if isreal(a), b = imag(b); end
if do_trans, b = b.'; end


