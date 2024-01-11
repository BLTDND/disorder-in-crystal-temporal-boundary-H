function [out] = norm_matrix(in)
    [r,c]=size(in);
    sum_in=0;
    for n=1:c
        sum_in=sum_in+abs(in(n))^2;
    end
    if sum_in == 0
        out = in;
    else
        out = in./sum_in;
    end
end