clear; close all; clc
fprintf('Loading and Visualizing Sources ... \n')
load('icaTest.mat'); % Original signals stored in matrices U (3x40) and A(3x3)
load('sounds.mat'); 
srcMat=[sounds(2,:); sounds(3,:); sounds(5,:)];
plot(0,0);
hold on;
xlabel('Time');
%ylabel('Frequency');
offSet=-1;
label='src';
offSet = addtoPlot(srcMat, offSet, label); 
fprintf('\n Showing Original Source. Press enter to continue.\n');
numSrc = size(srcMat,1);
A=rand(numSrc, numSrc);
X = A*srcMat;
label='mix';
offSet = addtoPlot(X, offSet, label); 
fprintf('\n Showing Mixed Signals. Press enter to continue.\n');
eta = 0.01;
eta0 = eta;
T=1000;
num_iter=10000;
W = rand(size(A)) ./ 10;
for i=0:num_iter,
	Y = W*X;			% predict source matrix based on guessed mix matrix
	[delW, mygrad] = gradient(eta, Y, W);	% gradient descent - shift by delta
	W = W + delW;			
	eta = eta0 / (1 + (i/T));	
end;
Y = W*X;
Y = (Y - min(min(Y))) ./ (max(max(Y)) - min(min(Y)));
Y2 = Y .* 2.0;
label='rec';
offSet = addtoPlot(Y2, offSet, label); 
print -dpng 'icaCheckImage.png';
hold off;
corrMat = correlations(srcMat,Y);
checkfile='icaCorrelations.txt';
printCorrs(corrMat, checkfile);
fprintf('\n Recovered Source Signals. Press enter to continue.\n');
pause;
clear; close all; clc
fprintf('Loading and Visualizing Sources ... \n')
load('icaTest.mat'); % Original signals stored in matrices U (3x40) and A(3x3)
load('sounds.mat'); % Original signals stored in sounds
srcMat=U;
plot(0,0);
hold on;
xlabel('Time');
%ylabel('Frequency');
offSet=-1;
label='src';
offSet = addtoPlot(srcMat, offSet, label); 
fprintf('\n Showing Original Source. Press enter to continue.\n');
numSrc = size(srcMat,1);
A=rand(numSrc, numSrc);
X = A*srcMat;
label='mix';
offSet = addtoPlot(X, offSet, label); 
fprintf('\n Showing Mixed Signals. Press enter to continue.\n');
eta = 0.01;
num_iter=100000;
W = rand(size(A));
for i=0:num_iter,
	Y = W*X;			
	delW = gradient(eta, Y, W);	
	W = W + delW;		
end;
Y = W*X;				
label='rec';
offSet = addtoPlot(Y, offSet, label); 
print -dpng "recoveredImage.png";
hold off;
fprintf('\n Recovered Source Signals. Press enter to continue.\n');
pause;
 function [grad, mygrad] = gradient(eta, Y, W)
delW = zeros(size(W));
Z = sigmoid(Y);
Id = eye(size(Y,1));
grad = eta * (Id + (1-2*Z)*Y') * W;
mygrad = eta * ((1-2*Z)*Y' + (size(Y,2)* pinv(W)') );
end
function [corrMat] = correlations(original, myrec)
numSrc = size(original,1);
corrMat = zeros(numSrc,numSrc);
for i=1:numSrc,
	for j=1:numSrc,
		corrMat(i,j) = corr2(original(i,:),myrec(j,:));
	end;
end;
end
function g= sigmoid(B, y)
g = zeros(size(y));
g = 1 ./ (1 + exp(- (B .* y)));
end
function [grad, delmyW, delb] = wgradientbeta(eta, kappa, b, Y, W)
delW = zeros(size(W));
B = repmat(b, 1, size(Y,2));
Z = sigmoidb(B,Y);
Id = eye(size(Y,1));
grad = eta * (Id + (B .* (1-2*Z))*Y') * W;
delmyW = eta * (((B .* (1-2*Z))*Y') +pinv(W)') ;
delb = kappa * (sum((Y .* (1-2*Z) *Y')'))' ;
end
