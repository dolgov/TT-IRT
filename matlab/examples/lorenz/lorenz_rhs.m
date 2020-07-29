% A function for computing the right-hand side of the Lorenz ODE
function [f] = lorenz_rhs(t,x,d)
x = reshape(x,[],d);
f = x;
for i=1:d
    if (i<d)
        xp1 = x(:, i+1);
    else
        xp1 = x(:, 1);
    end
    if (i>1)
        xm1 = x(:, i-1);
    else
        xm1 = x(:, d);
    end
    if (i>2)
        xm2 = x(:, i-2);
    else
        xm2 = x(:, d-1);
    end
    f(:,i) = (xp1 - xm2).*xm1 - x(:,i) + 8;
end
f = f(:);
end
