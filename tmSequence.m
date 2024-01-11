function tm_seq = tmSequence(n, a, b)
    tm_seq = zeros(1, n);
    tm_seq(1) = 0;
    tm_seq(2) = 1;
    
    for i = 3:n
        tm_seq(i) = mod(tm_seq(i-2) + tm_seq(i-1), 2);
    end
    for i=1:n
        if tm_seq(i)==0
            tm_seq(i)=a;
        else
            tm_seq(i)=b;
        end
    end
end