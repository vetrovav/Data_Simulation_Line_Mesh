
% This example use of the LineMesh function reads in a set of 33 bivariate
% flood discharge and flood volume data and uses 10 million random
% simulations the data-generated line mesh distribution with respect to
% estimating the joint non-exceedance probability with respect to the
% two bivariate x-y points 1.2,5.7 and 3.3, 9.2, with lambda set at 1.0.

% The two bivariate joint non-exceedance probabilities are obtained as
% 0.206 and 0.995, respectively.

% See the Excel spreadsheet Yue_input.xlsx for more data information.



Data=xlsread('Yue_input.xlsx'); 

lambda=1;

C=[1.2, 5.7; 3.3, 9.2];


K=10000000; % number of simulations




Upper=0;


% pp returns the p values (non-exceedance probabilities in this example
% because Upper is set to 0.

pp= LineMesh(Data,lambda,K,C,Upper)


        
        