function u1 = shock_filter(u, dt, step, u_origin)
%实现shock_filter求解

minmod = @(a,b) sign(max(a.*b,0)).*sign(a).*min(abs(a),abs(b));

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
 
    ux = (ue - uw)./2; uy = (un - us)./2;
    uxx = ue + uw - 2*u; uyy = un + us - 2*u;
    uxy = (une - unw - use + usw)./4;
    L = uw+ue+us+un-4*u;
   % L = (ux.^2.*uxx+2*ux.*uy.*uxy+uy.^2.*uyy)./(ux.^2+uy.^2);
    u = u - dt*sqrt(minmod(ue-u,u-uw).^2+minmod(un-u,u-us).^2).*L;
    err(i+1) = norm(u - u_origin)/norm(u_origin);

end

u1 = u;
figure 
plot(0:dt:step*dt,err);
title('shock filter');
xlabel('t');
ylabel('L2 error');

end
