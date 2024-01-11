function fib = fibonacci(n)
    if n <= 0
        error('Input must be a positive integer.');
    end
    
    fib = zeros(1, n);
    fib(1) = 1;
    if n > 1
        fib(2) = 1;
        for i = 3:n
            fib(i) = fib(i-1) + fib(i-2);
        end
    end
end