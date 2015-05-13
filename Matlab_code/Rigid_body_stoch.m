%SO(3) Stochastic "Lie" variational integrator based on "A Lie Group Variational Integrator for the Attitude Dynamics of a
%Rigid Body with Applications to the 3D Pendulum"

clear; close all; 

J = eye(3); J(1,1) = 3; J(2,2) = 2; J(3,3) = 1.5;                           %Inertia tensor
Jinv = inv(J);                                                              %Inverse of inertia tensor
Omega0 = [0.5;-0.5; 0.4];                                                   %Initial \Omega
R0 = eye(3);                                                                %Initial R
R = R0;
h = 0.01;                                                                   %Time step
T_end = 100;
t = 0:h:T_end;

N = length(t);

Pi_k = J*Omega0;
Pi_t = zeros(3, N);                                                         %Array to store the values of \Pi(t_k)
Fk = zeros(3,3);                                                            %F_k = R_k^T \dot{R}_k
M = [0;0;0];                                                                %M is a force 

sigma1 = [0.8; 0.2; 0.2];
sigma = sigma1*sigma1';                                                     %Variance
dW = (1/sqrt(h))*sigma*randn(3,1);                                          %Wiener process
err = zeros(2, N);                                                          %Array to record errors at each time step
M = cross(dW, Pi_k);
tic 
 
                               

for i=1:N
   
     M_old = M;
     Fk = RodSolve(h,J, Pi_k, M);                                           %Solve implicit equation using Rodrigues' formula(eqn 25)
     R = R*Fk;                                                              %R_{k+1} = R_k F_k
     dW = (1/sqrt(h))*sigma*randn(3,1);  
     M = cross(dW, Pi_k);                                                   %Update the force M_{k+1}
     Pi_k = (Fk')*Pi_k + (h/2)*(Fk')*M_old + (h/2)*M;                       %(eqn 24) \Pi_{k+1} = F_k^T \Pi_k + \frac{h}{2} F_k^T M_k + \frac{h}{2}M_{k+1}

     Pi_t(:,i) = Pi_k;
     Omega = Jinv*Pi_k;
     err(2, i) = norm( Pi_k'*cross(Omega, Pi_k));
     err(1,i)  = norm(eye(3) - R*R');                                       %\| I - RR^T \|_2
     
end

toc

figure
plot3(Pi_t(1,:),Pi_t(2,:),Pi_t(3,:))
title('Solution of Stochastic SO(3) Rigid Body')

figure 
plot(t, Pi_t(1,:), t, Pi_t(2,:), t, Pi_t(3,:))
title('Components of Angular Momentum')
legend('\Pi_1','\Pi_2', '\Pi_3')
xlabel('Time')

figure
plot(t, err(2,:))
title('Lie Group Integrator Error | I - R*R^T |')
xlabel('Time')
