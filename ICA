% Machine Learning HomeWork 2 - Independent Component Analysis
% This puts together the ica algorithm

%%Init
clear; close all; clc

%%=========Load and Sample Visualization

% Load Matrices U/srcMat and A; U is the original sound source, A is the mixing matrix
 
fprintf('Loading and Visualizing Sources ... \n')
load('icaTest.mat'); % Original signals stored in matrices U (3x40) and A(3x3)
load('sounds.mat'); % Original signals stored in sounds

%srcMat=U;
%srcMat=sounds;  	% Takes too much time to run all 5
%srcMat=sounds(1:3,:);	% Take only 3 signals
%srcMat=sounds(2:4,:);	% Take only 3 signals
%srcMat=sounds(3:5,:);	% Take only 3 signals
srcMat=[sounds(2,:); sounds(3,:); sounds(5,:)];

plot(0,0);
hold on;
xlabel('Time');
%ylabel('Frequency');
offSet=-1;
label='src';
offSet = addtoPlot(srcMat, offSet, label); 

fprintf('\n Showing Original Source. Press enter to continue.\n');
%pause;

%Mix and get X - mixed signals

numSrc = size(srcMat,1);
A=rand(numSrc, numSrc);
%fprintf('\n Actual A.\n');
%A
X = A*srcMat;
label='mix';
offSet = addtoPlot(X, offSet, label); 

fprintf('\n Showing Mixed Signals. Press enter to continue.\n');
%pause;

eta = 0.01;
eta0 = eta;
T=1000;
num_iter=10000;

%Make some random guess of mix-matrix inverse
W = rand(size(A)) ./ 10;

for i=0:num_iter,
	Y = W*X;			% predict source matrix based on guessed mix matrix
	[delW, mygrad] = gradient(eta, Y, W);	% gradient descent - shift by delta
	W = W + delW;			% update W
	%W = W + (mygrad * 0.001);			% update W
	eta = eta0 / (1 + (i/T));	% annealing - learning rate
	%if(mod(i,1000)==0),
	%	fprintf('runs %d \n',i);
	%	W
	%	corrMat = correlations(srcMat,Y)
	%	%fflush(stdout);
	%end;
end;

Y = W*X;				% predict source matrix based on guessed mix matrix

Y = (Y - min(min(Y))) ./ (max(max(Y)) - min(min(Y)));
%Y2=Y;
Y2 = Y .* 2.0;

% Display the obtained signals
%hold off;

%plot(0,0);
%xlabel('Time');
label='rec';
offSet = addtoPlot(Y2, offSet, label); 
%for i=1:numSrc,
%	plot(Y(i,:), sprintf('%s',plotColors(i)));
%	%print -dpng "recoveredImage1.png";
%	%fprintf('Check Plot. Press Enter to continue');
%	pause;
%end;
print -dpng 'icaCheckImage.png';
hold off;

% Compute correlation matrix to see which signals match and how well.
corrMat = correlations(srcMat,Y);
checkfile='icaCorrelations.txt';
printCorrs(corrMat, checkfile);

%save('235recovered.mat', 'Y2', 'corrMat');
%Check sounds
%soundsc(Y2(1,:),11025);
fprintf('\n Recovered Source Signals. Press enter to continue.\n');
pause;


%% Run choosing all sources in a loop

%for s1=1:3,
%	for s2=s1+1:4,
%		for s3=s2+1:5,
%			fprintf('Running on sources %d %d %d', s1,s2,s3);
%			srcMat=[sounds(s1,:); sounds(s2,:); sounds(s3,:)];
%			filename = sprintf('t-4kAnnealrecover-srcs-%d%d%d',s1,s2,s3);
%			plot(0,0);
%			hold on;
%			xlabel('Time');
%			offSet=-1;
%			title(sprintf('Graphs (Sources %d %d %d)',s1,s2, s3));
%			label='src';
%			offSet = addtoPlot(srcMat, offSet, label); 
%			numSrc = size(srcMat,1);
%			A=rand(numSrc, numSrc);
%			X = A*srcMat;
%			label='mix';
%			offSet = addtoPlot(X, offSet, label); 
%			eta = 0.01;
%			T=1000;
%num_iter=500;
%			eta0 = eta;
%			num_iter=100000;
%			W = rand(size(A)) ./ 10;
%			for i=0:num_iter,
%				Y = W*X;			
%				delW = gradient(eta, Y, W);	
%				W = W + delW;			
%				eta = eta0 / (1 + (i/T));	% annealing - learning rate
%			end;
%			Y = W*X;				
%			Y = (Y - min(min(Y))) ./ (max(max(Y)) - min(min(Y)));
%			Y2 = Y .* 2.0;
%			label='rec';
%			offSet = addtoPlot(Y2, offSet, label); 
%			imagefile = [filename, '.png'];
%			print('-dpng',imagefile);
%			hold off;
%			corrMat = correlations(srcMat,Y);
%			textfile = [filename, '.txt'];
%			printCorrs(corrMat, textfile);
%			matrixname = [filename, '.mat'];
%			save(matrixname, 'Y2', 'corrMat');
%		end;
%	end;
%end;
		
% Machine Learning HomeWork 2 - Independent Component Analysis
% This puts together the ica algorithm

%%Init
clear; close all; clc

%%=========Load and Sample Visualization

% Load Matrices U/srcMat and A; U is the original sound source, A is the mixing matrix
 
fprintf('Loading and Visualizing Sources ... \n')
load('icaTest.mat'); % Original signals stored in matrices U (3x40) and A(3x3)
load('sounds.mat'); % Original signals stored in sounds

srcMat=U;
%srcMat=sounds;  	% Takes too much time to run all 5
%srcMat=sounds(1:2,:);	% Take only 2 signals

plot(0,0);
hold on;
xlabel('Time');
%ylabel('Frequency');
offSet=-1;
label='src';
offSet = addtoPlot(srcMat, offSet, label); 

fprintf('\n Showing Original Source. Press enter to continue.\n');
%pause;

%Mix and get X - mixed signals

numSrc = size(srcMat,1);
A=rand(numSrc, numSrc);
%fprintf('\n Actual A.\n');
%A
X = A*srcMat;
label='mix';
offSet = addtoPlot(X, offSet, label); 

fprintf('\n Showing Mixed Signals. Press enter to continue.\n');
%pause;

eta = 0.01;
num_iter=100000;

%Make some random guess of mix-matrix inverse
W = rand(size(A));

for i=0:num_iter,
	Y = W*X;			% predict source matrix based on guessed mix matrix
	delW = gradient(eta, Y, W);	% gradient descent - shift by delta
	W = W + delW;			% update W
	%if(mod(i,1000)==0),
	%	fprintf('runs %d \n',i);
	%	corrMat = correlations(srcMat,Y);
	%	W
	%	fflush(stdout);
	%end;
end;

Y = W*X;				% predict source matrix based on guessed mix matrix

% Display the obtained signals

label='rec';
offSet = addtoPlot(Y, offSet, label); 
print -dpng "recoveredImage.png";
hold off;

% Compute correlation matrix to see which signals match and how well.
%corrMat = correlations(srcMat,Y);
%printCorrs(corrMat);

fprintf('\n Recovered Source Signals. Press enter to continue.\n');
pause;

 function [grad, mygrad] = gradient(eta, Y, W)
%Gradient function, it gives the value by which W must be updated to get new W
% \del W = eta (I + (1-2Z)Y')W where
%	eta - learning rate
%	Z = sigmoid(Y)
%	W is the predicted matrix

delW = zeros(size(W));

Z = sigmoid(Y);
Id = eye(size(Y,1));
grad = eta * (Id + (1-2*Z)*Y') * W;
mygrad = eta * ((1-2*Z)*Y' + (size(Y,2)* pinv(W)') );

end


function [corrMat] = correlations(original, myrec)
%Compute correlation between original signal and recovered signals
% Output matrix of correlation co-efficient of every pair of original and recovered signals.

numSrc = size(original,1);
corrMat = zeros(numSrc,numSrc);

for i=1:numSrc,
	for j=1:numSrc,
		corrMat(i,j) = corr2(original(i,:),myrec(j,:));
	end;
end;

end


function g= sigmoid(B, y)
% The sigmoid/logistic function
% g(y) = 1/(1 + e^{-by})
% Returns the sigmoid of y, where y can be a single value, vector or a matrix
% If b and y are single values we directly multiply them. If b is a scalar and y is a vector, we multiply the scalar with the vector.
% B is the same dimension as y!

g = zeros(size(y));
g = 1 ./ (1 + exp(- (B .* y)));

end

function [grad, delmyW, delb] = wgradientbeta(eta, kappa, b, Y, W)
%Gradient function, it gives the value by which W must be updated to get new W
% \del W = eta (I + (1-2Z)Y')W where
%	eta - learning rate
%	Z = sigmoid(Y)
%	W is the predicted matrix
%	b is a column vector with one value for each source

delW = zeros(size(W));

B = repmat(b, 1, size(Y,2));
Z = sigmoidb(B,Y);
Id = eye(size(Y,1));
grad = eta * (Id + (B .* (1-2*Z))*Y') * W;
delmyW = eta * (((B .* (1-2*Z))*Y') +pinv(W)') ;
delb = kappa * (sum((Y .* (1-2*Z) *Y')'))' ;

end


   
    
