% This is an implementation of Jacobi that can be run in Matlab's debugger.
% It is useful for debugging the C code.

% The problem
n = 8;
A = diag(2.^(-4:3));
y = ones(n, 1);
x = zeros(n, 1);

% The algorithm
newcoord = @(i, x)(x(i) + 1/A(i, i) * (y(i) - A(i, :) * x));
step = @(i, x)([x(1:i - 1); newcoord(i, x); x(i + 1:n)]);

for row = [7 7 7 7 3 7 4 4 7 1 7 2 7 7 3 7] % zero-based row indices
  x = step(row + 1, x);
  disp([sprintf('%5d:', row), sprintf('%6g', x), sprintf('; res = %6g', norm(y - A * x))]);
end
