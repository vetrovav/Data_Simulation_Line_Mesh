function  p = LineMesh(Data,lambda,K,C,Upper)
% Written by Earl Bardsley, University of Waikato
%
% This function uses K simulations from a line mesh distribution to estimate (as binomial proportions) either joint exceedance
% probabilities or joint non-exceedance probabilities, with respect to a set of N-variate specified points
% In this instance a finite mixture distribution of two normal distribution is applied for each line in the line mesh distribution
% (see paper for further description)
%
 

% Input
% Data: J x N array containing J recorded multivariate data values of dimension N
% lambda: smoothing parameter where 0.0001 <= lambda <= 0.9999
% K: user-specified number of simulations to be carried out
% C: R x N array containing R specified N-dimensional values for exceedance or non-exceedance probabilities
% Upper: flag set to zero or 1
%       Upper = 1 will cause the function to return a set of p values as joint
%       exceedance probabilities with respect to the R specified N-variate points in C
%       Upper = 0 will cause the function to return a set of p values as joint
%       non-exceedance probabilities with respect to the R specified N-variate points in C

% 		If Upper=1, the returned set of R values of p are estimates of the joint exceedance probabilities:
%       pr{X1 > C(1,1), X2 > C(1,2) ...XN > C(1,N)}
%		pr{X1 > C(2,1), X2 > C(2,2) ...XN > C(2,N)}
%       .....
%    	pr{X1 > C(R,1), X2 > C(R,2) ...XN > C(R,N)}

% 		If Upper=0, the returned set of R values of p are estimates of the joint non-exceedance probabilities:
%       pr{X1 < C(1,1), X2 < C(1,2) ...XN < C(1,N)}
%		pr{X1 < C(2,1), X2 < C(2,2) ...XN < C(2,N)}
%       .....
%    	pr{X1 < C(R,1), X2 < C(R,2) ...XN < C(R,N)}



D=size(Data);
J=D(1);
N=D(2);

BinCo=nchoosek(J,2);
v=1:J;
s=nchoosek(v,2);

% Select a random sample of K lines (which may be repeats) from the BinCo lines defined by all possible data point pairs
MM=randi(BinCo, K, 1);
%DataPoint1 is an array of size K,N
%DataPoint2 is an array of size K,N
DataPoint1=s(MM,1);
DataPoint2=s(MM,2);

DataValue1=Data(DataPoint1,:);
DataValue2=Data(DataPoint2,:);


% Find the maximum absolute difference between the respective coordinate pairs(to avoid small variances in the normal distributions)
difference=zeros(K,N);
n=1:N;
difference(:,n)=abs(DataValue1(:,n)-DataValue2(:,n));
[~, Ival]=max(difference,[],2); % Ival is a K*1 column vector containing the index values of the largest absolute difference


% For each sampled line, find the common standard deviation of the two normal distributions concerned, given the lambda value
% (see paper)

x=zeros(K,2);
for I=1:K
    x(I,1)=DataValue1(I,Ival(I));
    x(I,2)=DataValue2(I,Ival(I));
end

vmax=zeros(K,1); % the array vmax holds the maximum possible s^2 values
s=zeros(K,1);
vmax(:)=0.25.*(x(:,1).^2 + x(:,2).^2)-0.5.*(x(:,1).* x(:,2));
if lambda < 0.0001
	lambda=0.0001;
end
s(:)=sqrt(lambda.*vmax(:)); % normal distribution standard deviation values


% Find the corresponding two normal distribution mu values (m1 and mu2), given the lambda value

av=zeros(K,1);
mu1=zeros(K,1);
mu2=zeros(K,1);
delta=zeros(K,1);

av(:)=(x(:,1)+x(:,2))./2;

if lambda>0.9999 % set m1 = mu2 = av, if lambda is very near 1
      mu1(:)=av(:);
      mu2(:)=av(:);
else
      delta(:)=sqrt(0.5*(x(:,1).^2+x(:,2).^2)-av(:).^2-s(:).^2); 
      mu1(:)=av(:)-delta(:);
      mu2(:)=av(:)+delta(:);
end


% Select either mu1 or mu2 with equal probability 
rv=rand(K,1);
rvi=round(rv);
rvj=zeros(K,1);
rvj(:)=(rvi(:)+1).*(1-rvi(:)); 

mu=zeros(K,1);
mu(:)=mu1(:).*rvi(:)+mu2(:).*rvj(:);


% Generate K univariate normal random variables from the mu, s arrays
% normsim is an array of size K,1
normsim=normrnd(mu,s);



% Calculate the sets of "simulated" N values using linear interpolations,
% and store them in the array linesim
linesim=zeros(K,N);
yy1=zeros(K,1);
yy2=zeros(K,1);
yx=zeros(K,1);
xx1=x(:,1);
xx2=x(:,2);
yx(:)=normsim(:);

for I=1:N
	yy1(:)=DataValue1(:,I);
	yy2(:)=DataValue2(:,I);
      
	linesim(:,I)=yy1(:)+(yx(:,1)-xx1(:)).*(yy2(:)-yy1(:))./(xx2(:)-xx1(:));	
end

% Find the proportion of the simulated points in N-space which are greater than (or less than) the specified values.
sz=size(C);
R=sz(1);
p=zeros(R,1);


if Upper==1
	for I=1:R
		B = repmat(C(I,:),K,1);
		hh=linesim >B;
		kkount=min(hh,[],2);
		p(I)=sum(kkount);
	end
end

if Upper==0 
	for I=1:R
		B = repmat(C(I,:),K,1);
		hh=linesim <B;
		kkount=min(hh,[],2);
		p(I)=sum(kkount);
	end
end


	
p=p./K;	


end

    


