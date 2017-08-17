function LUTp11 = bivbern_lut(theta, lambda)
%BI_CALC_LUT Computes a look-up table for the joint firing probability

N1 = length(theta);
N2 = length(lambda);
% look-up table
LUTp11 = zeros(N1, N1, N2);
for k = 1:N2
    lam = lambda(k);
    for j = 1:N1
        th2 = theta(j);
        for i = j:N1
            p11 = bivnor(theta(i), th2, lam);
            LUTp11(i,j,k) = p11;
            if i ~= j
                LUTp11(j,i,k) = p11;
            end
        end
    end
end
end

