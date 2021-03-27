
function [Dk,i] = VAfind_Dk(f, tol, imax, D_init)


Dk_old = D_init ;
h = 1e-2;
Jacobi_old = ( f(Dk_old + ones(size(Dk_old)) * h )-f(Dk_old) )/h ;
Dk_new = Dk_old - Jacobi_old\f(Dk_old);

i = 0;

while (i<imax)
    if norm(f(Dk_new))<tol
        Dk = Dk_new;
        break
    else
        i = i+1;
        Dk_old = Dk_new; 
        Jacobi_old = (f(Dk_old+h)-f(Dk_old))/h;
        if cond(Jacobi_old) > 1e+16
            disp('Singular Jacobian in Newton method')
        end
        Dk_new = Dk_old - Jacobi_old\f(Dk_old);
    end
end

if (norm(f(Dk_new))>tol)
    disp('Qua so lan lap')
end

Dk = Dk_new ;

disp('Nghiem cua phuong trinh la: ')
Dk

disp('So lan lap la:')
i
