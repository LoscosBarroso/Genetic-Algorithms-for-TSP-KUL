function out = rand_int(m, n, range)
%  RAND_INT Generate matrix of uniformly distributed random integers.
%  OUT = RAND_INT(M,N) creates an M by N matrix of random integers 0 or 1
%  OUT = RAND_INT(M,N,RANGE) creates an M by N matrix of random integers.
%               RANGE can be a scalar or a 2-element vector. 
%                 * if RANGE is a scalar, random integers are generated in 
%                     the closed interval [0, RANGE-1]
%                 * if range is a 2-element vector, random integers are
%                      generated in the closed interval [RANGE(1), RANGE(2)]


if nargin == 2
    a = 0;
    b = 1;
else
    switch size(range, 2)
        case 1
            a = 0;
            b = range-1;
        case 2
            a = min(range);
            b = max(range);
        otherwise
            out = -1;
            return;
    end
end

out = zeros(m,n);
c = b-a+1;

out = floor(c .* rand(m,n)) + a;

end