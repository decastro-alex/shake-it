
%
% produces a realization of the process m(t) which satisfies:  
% 
%  dm = c^2 m dt + c \sum_{i=1}^3 \hat{e_i} m dW_i 
%

% physical parameters

mknot=[0 0 1];  % initial condition
c=0.1;   

% numerical parameters

T=1000;
dt=0.1;
num_steps=ceil(T/dt);

% preprocessing

e1hat=[0 0 0; 0 0 -1; 0 1 0];
e2hat=[0 0 1; 0 0 0; -1 0 0];
e3hat=[0 -1 0; 1 0 0; 0 0 0];

% allocate space

t_vec=zeros(num_steps,1);
m_vec=zeros(num_steps,3);
m_vec(1,:)=mknot;

for i=2:num_steps
    m0=m_vec(i-1,:)';
    m1=expm(c*sqrt(dt)*(e1hat*randn+e2hat*randn+e3hat*randn))*m0;
    m_vec(i,:)=m1';
    t_vec(i)=t_vec(i-1)+dt;    
end

%% graphical outputs

figure(1);
plot3(m_vec(:,1),m_vec(:,2),m_vec(:,3),'k','LineWidth',2);
axis(1.5*[-1 1 -1 1 -1 1]);
xlabel('$x$','fontsize',20,'Interpreter','latex');
ylabel('$y$','fontsize',20,'Interpreter','latex');
zlabel('$z$','fontsize',20,'Interpreter','latex');
grid on;
box on;

figure(2); hold on;
plot(t_vec,m_vec(:,1),'color',[0 0 0],'LineWidth',2);
plot(t_vec,m_vec(:,2),'color',[0.45 0.45 0.45],'LineWidth',2);
plot(t_vec,m_vec(:,3),'color',[0.75 0.75 0.75],'LineWidth',2);
xlabel('$t$','fontsize',20,'Interpreter','latex');
ylabel('$m_i(t)$','fontsize',20,'Interpreter','latex');
lh=legend({'$m_1(t)$', '$m_2(t)$', '$m_3(t)$'}, 'location', 'best', 'Interpreter','latex', 'fontsize',20);
set(lh,'fontsize',20);
grid on;
box on;

figure(3);
plot(t_vec,m_vec(:,1).^2+m_vec(:,2).^2+m_vec(:,3).^2,'k','LineWidth',2);
xlabel('$t$','fontsize',20,'Interpreter','latex');
ylabel('$|m(t)|^2$','fontsize',20,'Interpreter','latex');
ylim([1/2 3/2]);
grid on;
box on;