function u1 = heat(u, dt, step, u_origin, nu)
%实现热方程的求解
if nargin < 5
    nu = 1;
end
err = zeros(1,step+1);
err(1) = norm(u - u_origin)/norm(u_origin);
for i = 1:step
    uw = u; uw(:,1:size(u,2)-1)=u(:,2:size(u,2));
    ue = u; ue(:,2:size(u,2))=u(:,1:size(u,2)-1);
    un = u; un(1:size(u,1)-1,:)=u(2:size(u,1),:);
    us = u; us(2:size(u,1),:)=u(1:size(u,1)-1,:);

    u = u + nu*dt*(uw + ue + us + un - 4*u);
    err(i+1) = norm(u - u_origin)/norm(u_origin);
end 
u1 = u;
figure 
plot(0:dt:step*dt,err);
title('heat PDE');
xlabel('t');
ylabel('L2 error');

end
