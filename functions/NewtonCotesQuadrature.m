function I = NewtonCotesQuadrature(a,b,y)
%   Summary of this function goes here.
%   I = NewtonCotesQuadrature(a,b,y)
%   Detailed explanation goes here.
%   The function is for calculating the numerical integral for the equally
%   spaced discrete point by Newton-Cotes method.
%
%   IN      
%             a: the lower bound of the integral range. 
%             b: the upper bound of the integral range.
%             y: it's a vector, the function value at the discrete points. 
%
%  OUT   
%             I: the integral value.
%  
%  References:
%  You can find the proper formula from the topic of numerical computation.
%
%  Author(s): Yan Yingwei
%  Copyright: 2020-2025 
%  Revision: 1.0  Date: 9/27/2020
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology (SUSTech).

C = zeros(10,11);
C(1,1:2) = 1/2*[1 1];
C(2,1:3) = 1/6*[1 4 1];
C(3,1:4) = 1/8*[1 3 3 1];
C(4,1:5) = 1/90*[7 32 12 32 7];
C(5,1:6) = 1/228*[19 75 50 50 75 19];
C(6,1:7) = 1/840*[41 216 27 272 27 216 41];
C(7,1:8) = 1/17280*[751 3577 1323 2989 2989 1323 3577 751];
C(8,1:9) = 1/28350*[989 5888 -928 10496 -4540 10496 -928 5888 989];
C(9,1:10) = 1/89600*[2857 15741 1080 19344 5778 5778 19344 1080 15741 2857];
C(10,1:11) = 1/5996752*[16067 106300 -48525 272400 -260550 427368 -260550 272400 -48525 106300 16067];

l = length(y);
n = l-1;

I = 0;
for i=1:l
    I = I+C(n,i)*y(i);
end
I = (b-a)*I;
end

