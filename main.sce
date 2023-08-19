// code of zain ul abideen
clear; clc;

// setting up trivial values for constants.
m = 1; h = 1;
a = 10; n = 500; d = a ./ n;

// setting up the domain of our calculation.
x = [-a ./ 2 : d : +a ./ 2]

// defining a trivial harmonic oscillator potential.
function [y] = V(x)
    y = x .^ 2
endfunction

// constructing the first matrix. (didn't work)
//for i = 1 : n - 1
//    A(i, i) = -1 + ((V(x(i)) .* 2 .* m .* d .* d) ./ (h .* h))
//    
//    if (i + 1) <= (n - 1) then
//        A(i, i + 1) = 2
//    end
//    if (i + 2) <= (n - 1) then
//        A(i, i + 2) = -1
//    end
//end


// constructing the second matrix. (worked flawlesly)
for i = 1 : n - 1
    A(i, i) = 2 .* (1 + (m .* d .* d ./ (h .* h)) .* V(x(i)))
    
    if i < (n - 1) then
        A(i, i + 1) = -1
        A(i + 1, i) = -1
    end
end

// extracting the eigenvalues and eigenvectors.
[X, lambda] = spec(A)

// plotting the results.
plot(x', [0; X(:, 1); 0], '-r')
plot(x', [0; X(:, 2); 0], '--g')
plot(x', [0; X(:, 3); 0], ':') // and so on.
legend('ground state', 'first excited state', 'second excited state')

// see the relevant paper for more insight.

