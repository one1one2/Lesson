function u1 = pm(u, dt, step, u_origin);
%实现pm方程求解

c = @(s) 1./sqrt(s + 1);

err = zeros(1,step+1);
err(1) = norm(u - u_origin)/norm(u_origin);
for i = 1:step
    uw = u; uw(:,1:size(u,2)-1)=u(:,2:size(u,2));
    ue = u; ue(:,2:size(u,2))=u(:,1:size(u,2)-1);
    un = u; un(1:size(u,1)-1,:)=u(2:size(u,1),:);
    us = u; us(2:size(u,1),:)=u(1:size(u,1)-1,:);

    unw = un; unw(:,1:size(u,2)-1)=un(:,2:size(u,2));
    usw = us; usw(:,1:size(u,2)-1)=us(:,2:size(u,2));
    une = un; une(:,2:size(u,2))=un(:,1:size(u,2)-1);
    use = us; use(:,2:size(u,2))=us(:,1:size(u,2)-1);
    
    bp0 = c((ue-u).^2 + 1./16*(une+un-us-use).^2); 
    bm0 = c((uw-u).^2 + 1./16*(unw+un-us-usw).^2); 
    b0p = c((un-u).^2 + 1./16*(une+ue-uw-unw).^2); 
    b0m = c((us-u).^2 + 1./16*(use+ue-uw-usw).^2);

    u = u + dt*(bp0.*ue+bm0.*uw+b0p.*un+b0m.*us-(b0p+bp0+b0m+bm0).*u);   
    err(i+1) = norm(u - u_origin)/norm(u_origin);
end 
u1 = u;
figure 
plot(0:dt:step*dt,err);
title('Perona-Malik');
xlabel('t');
ylabel('L2 error');
end 
