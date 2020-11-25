function y = twopBVP(fvec, alpha, beta, L, N)

    dx = L/(N+1);
    
    %only considering the simplest problem: F(x,y) = f(x)
    %meaning that the right side is independet of y
    
    W = diag(-2*ones(1,N)) + diag(1*ones(1,N-1),1) + diag(1*ones(1,N-1),-1);
    W = W*(1/(dx)^2);
    
    fvec(1) = fvec(1) - alpha/(dx)^2;
    fvec(end) = fvec(end) - beta/(dx)^2;
       
    fvect = fvec';
    
    y = W\fvect;
    
end
