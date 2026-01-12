function [k1, k2, l_star] = kappa_search_doubly(d, rho)

    if d(end) / d(1) <= rho  
         k1 = 0; k2 = 0; 
         l_star = 0; 
         return;
    end

    n = length(d);

    head_sum = cumsum(d);
    tail_sum = cumsum(d, 'reverse'); 

    lowerbound = d(1); 

    k1 = 1;
    k2 = 2;

    while k2 <= n && d(k2)/rho <= d(k1)
        k2 = k2 + 1;
    end
    
    while k1 < n && k2 <= n
        
        if d(k2) / rho < d(k1 + 1)
            upperbound = d(k2) / rho;
            raise_k2 = true;
        else
            upperbound = d(k1 + 1);
            raise_k2 = false;
        end
        
        l_star = (head_sum(k1) + rho * tail_sum(k2)) / (k1 + rho^2 * (n - k2 + 1));
        
        if l_star > lowerbound && l_star <= upperbound 
            return;
        end
        
        lowerbound = upperbound;
        
        if raise_k2
            k2 = k2 + 1;
        else
            k1 = k1 + 1;
            
            while k2 <= n && d(k2)/rho <= d(k1)
                k2 = k2 + 1;
            end
        end
    end
end
