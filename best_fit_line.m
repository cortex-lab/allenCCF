function [m,p,s] = best_fit_line(x,y,z)

% x,y,z are n x 1 column vectors of the three coordinates
% of a set of n points in three dimensions. The best line,
% in the minimum mean square orthogonal distance sense,
% will pass through m and have direction cosines in p, so
% it can be expressed parametrically as x = m(1) + p(1)*t,
% y = m(2) + p(2)*t, and z = m(3)+p(3)*t, where t is the
% distance along the line from the mean point at m.
% s returns with the minimum mean square orthogonal
% distance to the line.
% RAS - March 14, 2005

[n,mx] = size(x); [ny,my] = size(y); [nz,mz] = size(z);
if (mx~=1)||(my~=1)||(mz~=1)||(ny~=n)||(nz~=n)
 error('The arguments must be column vectors of the same length.')
end

% m is the mean value of each dimension
m = [mean(x),mean(y),mean(z)];

% w is the mean-subtracted data
w = [x-m(1),y-m(2),z-m(3)]; % Use "mean" point as base

a = (1/n)*w'*w; % 'a' is a positive definite matrix of the data
[u,d,~] = svd(a); % 'eig' & 'svd' get same eigenvalues for this matrix
p = u(:,1)'; % Get eigenvector for largest eigenvalue
s = d(2,2)+d(3,3); % Sum the other two eigenvalues
