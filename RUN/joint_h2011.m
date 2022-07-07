function mutu = joint_h2011 (dif,spc,row,col)

  
    dif1 = dif;   

    spc1 = spc; 


%% To use in rounding after this
y = max(dif1);
z = max(spc1);


%% Change back to width x height

dif2 = reshape (dif1,row,col);
spc2 = reshape (spc1,row,col);


%% Rounding to 256 , grayscale image.
dif3 = round ( dif2*255/y);
spc3 = round ( spc2*255/z);


%% Build Joint Histogram of Diffuse and Specular
rows=size(dif3,1);
cols=size(dif3,2);
N=256;
h=zeros(N,N);
for m=1:rows;    %  col 
  for n=1:cols;  %  rows
    h(dif3(m,n)+1,spc3(m,n)+1)= h(dif3(m,n)+1,spc3(m,n)+1)+1;
  end
end


%% Built Histogram for Diffuse and Specular
dif4 = dif3(:)';
spc4 = spc3(:)';
c = 1:1:256 ;
dif5 = hist (dif4,c);
spc5 = hist (spc4,c);


%% Mutual Information,I using three parameters
mutu = 0;
I=0;
for p = 1:1:256
    for q = 1:1:256   
        if (h (p,q)==0 | dif5 (p)==0 | spc5 (q)==0)       
            I = 0;
        else
            dif6 (p) = dif5 (p) /  (rows*cols)  ;
            spc6 (q) = spc5 (q) /  (rows*cols)  ;
            h2 (p,q) = h (p,q)  /  (rows*cols)  ;
            I =  (h2 (p,q))* log ((h2 (p,q))/((dif6 (p))*(spc6 (q)))) ;
        end
     mutu = mutu + I;
    end   
end

%% Return MUTU