function x=tfsynthesis(timefreqmat,swin,timestep,numfreq)

%time-frequency synthesis
%TIMEFREQMAT is the complex matrix time-freq representation
%SWIN is the synthesis window
%TIMESTEP is the # of samples between adjacent time windows.
%NUMFREQ is the # of frequency components per time point.
%X contains the reconstructed signal.

swin=swin(:); %make synthesis window go column-wise
winlen=length(swin);
[numfreq numtime]=size(timefreqmat);
ind=rem((1:winlen)-1,numfreq)+1;
x=zeros((numtime-1)*timestep+winlen,1);

for i=1:numtime%overlap,window,andadd
    temp=numfreq*real(ifft(timefreqmat(:,i)));
    sind=((i-1)*timestep);
    rind=(sind+1):(sind+winlen);
    x(rind)=x(rind)+temp(ind).*swin;
end

function tfmat=tfanalysis(x,awin,timestep,numfreq)

%time-frequency analysis
%X is the time domain signal
%AWIN is an analysis window
%TIMESTEP is the # of samples between adjacent time windows.
%NUMFREQ is the # of frequency components per time point.
%TFMAT complex matrix time-freq representation

x=x(:); 
awin=awin(:); %make inputs go column-wise
nsamp=length(x);
wlen=length(awin);%calc size and init output t-f matrix
numtime=ceil((nsamp-wlen+1)/timestep);
tfmat=zeros(numfreq,numtime+1);

% so this loop samples x, I though x is already a sample.
for i=1:numtime
    sind=((i-1)*timestep)+1; % current start index
    tfmat(:,i)=fft(x(sind:(sind+wlen-1)).*awin,numfreq);
end
% below is inelegant... but apparently works.
i=i+1;
sind=((i-1)*timestep)+1;
lasts=min(sind,length(x));
laste=min((sind+wlen-1),length(x));
tfmat(:,end)=fft([x(lasts:laste);
                    zeros(wlen-(laste-lasts+1),1)].*awin,numfreq);
                    
                    
function smat=twoDsmooth(mat,ker)
%TWO2SMOOTH?Smooth 2D matrix.
%
%smat=twoDsmooth(mat,ker)
%
%MAT is the 2D matrix to be smoothed.
%KER is either
%(1)a scalar,in which case a ker?by?ker matrix of 1/ker?2 is used 
%   as the matrix averaging kernel
%(2)a matrix which is used as the averaging kernel.
%
%SMATisthesmoothedmatrix(samesizeasmat).
if numel(ker)==1,%prod(size(ker))==1 if ker is a scalar
    kmat=ones(ker,ker)/ker^2;
else
    kmat=ker;
end
%make kmat have odd dimensions
[kr kc]=size(kmat); 
if rem(kr,2)==0, %Remainder after division
    kmat=conv2(kmat,ones(2,1))/2;
    kr=kr+1;
end
if rem(kc,2)==0,
    kmat=conv2(kmat,ones(1,2))/2;
    kc=kc+1;
end
[mr mc]=size(mat);
fkr=floor(kr/2);%number of rows to copy on top and bottom
fkc=floor(kc/2);%number of columns to copy on either side
smat=conv2(...
    [mat(1,1)*ones(fkr,fkc) ones(fkr,1)*mat(1,:)...
    mat(1,mc)*ones(fkr,fkc);
    mat(:,1)*ones(1,fkc) mat mat(:,mc)*ones(1,fkc)
    mat(mr,1)*ones(fkr,fkc) ones(fkr,1)*mat(mr,:)...
    mat(mr,mc)*ones(fkr,fkc)],...
    rot90(kmat,2),'valid'); %flipud(fliplr(kmat)) replaced by rot90(kmat,2)

clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    setp 1,2,3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. analyze the signals - STFT
%1) Create the spectrogram of the Left and Right channels.
wlen=1024;
timestep=512;
numfreq=1024;
awin=hamming(wlen);%analysis window is a Hamming window Looks like Sine on [0,pi]
[x1,fs] = audioread('data/x1_reverb.wav');
x2 = audioread('data/x2_reverb.wav');
tf1=tfanalysis(x1,awin,timestep,numfreq);%time-freq domain
tf2=tfanalysis(x2,awin,timestep,numfreq);%time-freq domain

tf1(1,:)=[];
tf2(1,:)=[];%remove dc component from mixtures
%eps is the a small constant to avoid dividing by zero frequency in the delay estimation

%calculate pos/neg frequencies for later use in delay calc ??
freq=[(1:numfreq/2) ((-numfreq/2)+1:-1)]*(2*pi/(numfreq)); % freq looks like saw signal
fmat=freq(ones(size(tf1,2),1),:)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%2.calculate alpha and delta for each t-f point
%2) For each time/frequency compare the phase and amplitude of the left and
%   right channels. This gives two new coordinates, instead of time-frequency 
%   it is phase-amplitude differences.
R21=(tf2+eps)./(tf1+eps);%time-freqratioofthemixtures

%%%2.1HERE WE ESTIMATE THE RELATIVE ATTENUATION (alpha)%%%
a=abs(R21);%relative attenuation between the two mixtures
alpha=a-1./a;%'alpha' (symmetric attenuation)
%%%2.2HERE WE ESTIMATE THE RELATIVE DELAY (delta)%%%%
delta=-imag(log(R21))./fmat;% imaginary part, 'delta' relative delay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3.calculate weighted histogram
%3) Build a 2-d histogram (one dimension is phase, one is amplitude) where 
%   the height at any phase/amplitude is the count of time-frequency bins that
%   have approximately that phase/amplitude.
p=1; q=0; %powers used to weight histogram
tfweight=(abs(tf1).*abs(tf2)).^p.*abs(fmat).^q; %weights vector
maxa=0.7;
maxd=3.6;%histogram boundaries for alpha, delta

abins=35;
dbins=50;%number of hist bins for alpha, delta

%only consider time-freq points yielding estimates in bounds
amask=(abs(alpha)<maxa)&(abs(delta)<maxd);
alphavec=alpha(amask);
deltavec=delta(amask);
tfweight=tfweight(amask);

%determine histogram indices (sampled indices?)
alphaind=round(1+(abins-1)*(alphavec+maxa)/(2*maxa));
deltaind=round(1+(dbins-1)*(deltavec+maxd)/(2*maxd));

%FULL-SPARSE TRICK TO CREATE 2D WEIGHTED HISTOGRAM
%A(alphaind(k),deltaind(k)) = tfweight(k), S is abins-by-dbins
A=full(sparse(alphaind,deltaind,tfweight,abins,dbins));
%smooththehistogram-localaverage3-by-3neighboringbins
A=twoDsmooth(A+10,1);

%plot2-Dhistogram
mesh(linspace(-maxd,maxd,dbins),linspace(-maxa,maxa,abins),A);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    step 4,5,6,7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.peak centers (determined from histogram) THIS IS DONE BY HUMAN.
%4) Determine how many peaks there are in the histogram.
%5) Find the location of each peak. 

numsources=5;

% peak estimates provided in source code
peakdelta=[-2 -2 0 2 2];
peakalpha=[.19 -.21 0 .19 -.21];

% my own peak estimates. I'm not sure the ones provided are correct.
peakdelta=[-1.4 .66 1.25 1.9 .51];
peakalpha=[.4 -.29 .12 .66 -.5];

%convert alpha to a
peaka=(peakalpha+sqrt(peakalpha.^2+4))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5.determine masks for separation
%6) Assign each time-frequency frame to the nearest peak in phase/amplitude 
%  space. This partitions the spectrogram into sources (one peak per source)

bestsofar=Inf*ones(size(tf1));
bestind=zeros(size(tf1));
for i=1:length(peakalpha)
    score=abs(peaka(i)*exp(-sqrt(-1)*fmat*peakdelta(i))...
        .*tf1-tf2).^2/(1+peaka(i)^2);
    mask=(score<bestsofar);
    bestind(mask)=i;
    bestsofar(mask)=score(mask);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6.&7.demix with ML alignment and convert to time domain
%7) Then you create a binary mask (1 for each time-frequency point belonging to my source, 0 for all other points)
%8) Mask the spectrogram with the mask created in step 7.
%9) Rebuild the original wave file from 8.
%10) Listen to the result.
est=zeros(numsources,length(x1));%demixtures
for i=1:numsources
    mask=(bestind==i);
    esti=tfsynthesis([zeros(1,size(tf1,2));
        ((tf1+peaka(i)*exp(sqrt(-1)*fmat*peakdelta(i)).*tf2)...
            ./(1+peaka(i)^2)).*mask],...
             sqrt(2)*awin/1024,timestep,numfreq);
    est(i,:)=esti(1:length(x1))';

    %add back into the demix a little bit of the mixture
    %as that eliminates most of the masking artifacts
    soundsc(est(i,:)+0.05*x1',fs);% original code seems to have missed the transpose play demixture
end



function [k1,kn,err] = decompose_kernel(h_orig)
% This function does the decomposition of a separable nD kernel into
% its 1D components, such that a convolution with each of these
% components yields the same result as a convolution with the full nD
% kernel, at a drastic reduction in computational cost.
%
% SYNTAX:
% =======
%    [K1,KN,ERR] = DECOMPOSE_KERNEL(H)
% computes a set of 1D kernels K1{1}, K1{2}, ... K1{N} such that the
% convolution of an image with all of these in turn yields the same
% result as a convolution with the N-dimensional input kernel H:
%    RES1 = CONVN(IMG,H);
%    RES2 = IMG;
%    FOR II=1:LENGTH(K1)
%       RES2 = CONVN(RES2,K1{II});
%    END
%
% KN is the reconstruction of the original kernel H from the 1D
% kernels K1, and ERR is the sum of absolute differences between H
% and KN.
% The syntax mimics Dirk-Jan Kroon's submission to the FileExchange
% (see below).
%
% EXPLANATION:
% ============
%
% In general, for a 2D kernel H, the convolution with 2D image F:
%    G = F * H
% is identical to the convolution of the image with column vector H1
% and convolution of the result with row vector H2:
%    G = ( F * H1 ) * H2   .
% In MATLAB speak, this means that
% > CONV2(F,H) == CONV2(CONV2(F,H1),H2)
%
% Because of the properties of the convolution,
%    ( F * H1 ) * H2 = F * ( H1 * H2 )   ,
% meaning that the convolution of the two 1D filters with each other
% results in the original filter H. And because H1 is a column vector
% and H2 a row vector,
%    H = H1 * H2 = H1 H2   .
% Thus, we need to find two vectors whose product yields the matrix H.
% In MATLAB speak we need to solve the equation
% > H1*H2 == H
%
% The function in the standard MATLAB toolbox, FILTER2, does just
% this, and it does it using singular value decomposition:
%    U S V' = H   ,
%    H1(i) = U(i,1)  S(1,1)^0.5   ,
%    H2(i) = V(i,1)* S(1,1)^0.5   .  (the * here is the conjugate!)
%
% Note that, if the kernel H is separable, all values of S are zero
% except S(1,1). Also note that this is an under-determined problem,
% in the sense that
%    H = H1 H2 = ( a H1 ) ( 1/a H2 )   ;
% that is, it is possible to multiply one 1D kernel with any value
% and compensate by dividing the other kernel with the same value.
% Our solution will, in effect, just choose one of the infinite number
% of (equivalent) solutions.
%
% To extend this concept to nD, what we need to understand is that it
% is possible to collapse all dimensions except one, obtaining a 2D
% matrix, and solve the above equation. This results in a 1D kernel
% and an (n-1)D kernel. Repeat the process until all you have is a
% set of 1D kernels and you're done!
%
% This function is inspired by a solution to this problem that
% Dirk-Jan Kroon posted on the File Exchange recently:
% http://www.mathworks.com/matlabcentral/fileexchange/28218-separate-kernel-in-1d-kernels
% His solution does the whole decomposition in one go, by setting up
% one big set of equations. He noted a problem with negative values,
% which produce complex 1D kernels. The magnitude of the result is
% correct, but the sign is lost. He needs to resort to some heuristic
% to determine the sign of each element. What he didn't notice (or
% didn't mention) is the problem that his solution has with 0 values.
% The SVD solution doesn't have this problem, although it sometimes
% does produce a slightly worse solution. For example, in the first
% example below, Dirk-Jan Kroon's solution is exact, whereas this one
% produces a very small error. Where Dirk-Jan Kroon's solution cannot
% find the exact solution, this algorithm generally does better.
% 
% EXAMPLES:
% =========
%
% Simplest 5D example:
%
%    H = ones(5,7,4,1,5);
%
%    [K1,~,err] = SeparateKernel(H); % D.Kroon's submission to FileEx.
%    err
%
%    [k1,~,err] = decompose_kernel(H);
%    err
%
% 2D example taken from Dirk-Jan Kroon's submission:
%
%    a = permute(rand(4,1),[1 2 3])-0.5;
%    b = permute(rand(4,1),[2 1 3])-0.5;
%    H = repmat(a,[1 4]).*repmat(b,[4 1]);
%
%    [K1,~,err] = SeparateKernel(H);
%    err
%
%    [k1,~,err] = decompose_kernel(H);
%    err
%
% 2D example for which Dirk-Jan Kroon's solution has problems:
%
%    H = [1,2,3,2,1]'*[1,1,3,0,3,1,1];
%   
%    [K1,~,err] = SeparateKernel(H);
%    err
%   
%    [k1,~,err] = decompose_kernel(H);
%    err
%
% 3D example that's not separable:
%
%    H = rand(5,5,3);
%
%    [K1,~,err] = SeparateKernel(H);
%    err
%
%    [k1,~,err] = decompose_kernel(H);
%
% Example to apply a convolution using the decomposed kernel:
%
%    img = randn(50,50,50);
%    h = ones(7,7,7);
%    tic;
%    res1 = convn(img,h);
%    toc
%    k1 = decompose_kernel(h);
%    tic;
%    res2 = img;
%    for ii=1:length(k1)
%       res2 = convn(res2,k1{ii});
%    end
%    toc
%    rms_diff = sqrt(mean((res1(:)-res2(:)).^2))

% Copyright 2010, Cris Luengo,
%   Centre for Image Analysis,
%   Swedish University of Agricultural Sciences and Uppsala University,
%   Uppsala, Sweden.
% 20 July 2010.

% Save the input kernel to compare to later on
h = h_orig;        
% Create cell array for output
n = ndims(h);
k1 = cell(1,n);
% Decompose last dimension iteratively until we've got a 2D kernel left
s = size(h);
for ii=n:-1:3
   h = reshape(h,[],s(end));       % collapse all dims except last one
   [h,k1{ii}] = decompose_2D(h);
   s = s(1:end-1);
   h = reshape(h,s);               % restore original shape with 1 fewer dim
   k1{ii} = shiftdim(k1{ii},2-ii); % the 1D kernel must be in the right shape
end
% Decompose the final 2D kernel
[k1{1},k1{2}] = decompose_2D(h);


% Reconstruct nD matrix from 1D components
kn = k1{1}*k1{2};
for ii=3:length(k1)
   s = ones(ii,1);
   s(ii) = length(k1{ii});
   kn = repmat(kn,s);
   s = size(kn);
   s(end+1:ii) = 1; % in case kn has ending singleton dimensions, which will be ignored
   s(ii) = 1;
   kn = kn.*repmat(k1{ii},s);
end

% Calculate the error we made
err = sum(abs(h_orig(:)-kn(:)));


function [h1,h2] = decompose_2D(h)
% This does the decomposition of a 2D kernel into 2 1D kernels
% More or less a direct copy-paste from the FILTER2 function
% in the standard MATLAB toolbox.
[ms,ns] = size(h);
if (ms == 1)
   h1 = 1;
   h2 = h;
elseif (ns == 1)
   h1 = h;
   h2 = 1;
else
   separable = false;
   if all(isfinite(h(:)))
      % Check rank (separability) of kernel
      [u,s,v] = svd(h);
      s = diag(s);
      tol = length(h) * eps(max(s));
      rank = sum(s > tol);   
      separable = rank==1;
   end
   if ~separable
      error('Sorry, kernel is not separable.');
   end
   h1 = u(:,1) * sqrt(s(1));
   h2 = (conj(v(:,1)) * sqrt(s(1)))';
end

