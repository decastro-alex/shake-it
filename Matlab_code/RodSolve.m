function F = RodSolve(h, J, Pik, M)
    
    F1 = [0.1;0.1;0.1];
    
%     F_n = norm(F1);
%     
%     
%     Fk = (sin(F_n)/F_n)*(J*F1) + ((1 - cos(F_n))/(F_n^2))*(cross(F1, J*F1)) - (h*Pik + (h^2/2)*M);
%     
    count = 1;
    tol = 0.00001; err =1;
    MAX_ITER = 100;
    
    while((count < MAX_ITER)&&(err > tol))
        F_n = norm(F1);
        Fk = (sin(F_n)/F_n)*(J*F1) + ((1 - cos(F_n))/(F_n^2))*(cross(F1, J*F1)) - (h*Pik + (h^2/2)*M);
        Fd = (sin(F_n)/F_n)*J + ((cos(F_n)*F_n - sin(F_n))/(F_n^3))*(J*F1*F1') + ((sin(F_n)*F_n -2*(1 - cos(F_n)))/(F_n^4))*((cross(F1, J*F1)*F1')) + ((1 - cos(F_n))/(F_n^2))*(-skew(-J*F1) + skew(F1)*J);
        V = -inv(Fd)*Fk;
        
        F1 = F1 + V;
        count = count + 1;
        err = norm(V);
    end
%     count
%     err
    SFk = skew(F1);
    
    F = eye(3) + (sin(F_n)/F_n)*SFk + ((1 - cos(F_n))/(F_n^2))*(SFk*SFk);
end