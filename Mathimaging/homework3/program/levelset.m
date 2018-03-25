function u1 = levelset(I, dt, step, alpha, tol, Reinit_step);
%实现pm方程求解

scale = 1e5;

fun_g = @(s) 1./(1 + scale*s.^2);

m = size(I,1); 
n = size(I,2);

Iw = I; Iw(:,1:n-1)=I(:,2:n);
Ie = I; Ie(:,2:n)=I(:,1:n-1);
In = I; In(1:m-1,:)=I(2:m,:);
Is = I; Is(2:m,:)=I(1:m-1,:);
    
Ix = (Ie - Iw)/2;
Iy = (In - Is)/2;
dI = (Ix.^2 + Iy.^2).^0.5;
    
g = fun_g(dI);

gw = g; gw(:,1:n-1)=g(:,2:n);
ge = g; ge(:,2:n)=g(:,1:n-1);
gn = g; gn(1:m-1,:)=g(2:m,:);
gs = g; gs(2:m,:)=g(1:m-1,:);
    
gx = (ge - gw)/2;
gy = (gn - gs)/2;

u = zeros(size(I));
for i = 1:m
  for j = 1:n
    u(i,j) = -min(min(abs(i-1),abs(m-i)), min(abs(j-1),abs(n-j)));
  end
end

for i = 1:step
    uw = u; uw(:,1:n-1)=u(:,2:n);
    ue = u; ue(:,2:n)=u(:,1:n-1);
    un = u; un(1:m-1,:)=u(2:m,:);
    us = u; us(2:m,:)=u(1:m-1,:);

    unw = un; unw(:,1:n-1)=un(:,2:n);
    usw = us; usw(:,1:n-1)=us(:,2:n);
    une = un; une(:,2:n)=un(:,1:n-1);
    use = us; use(:,2:n)=us(:,1:n-1);
  
    ux = (ue - uw)/2;
    uy = (un - us)/2;
    du = sqrt(ux.^2 + uy.^2);
    uxx = (ue + uw - 2*u);
    uyy = (un + us - 2*u);
    uxy = (une - unw - use + usw)/4;
    % \nabla^+u , \nabla^-u
    uxp = ue - u; uxm = u - uw;
    uyp = un - u; uym = u - us;
    dpu = (max(uxm,0).^2 + min(uxp,0).^2 + max(uym,0).^2 + min(uyp,0).^2).^0.5;
    dmu = (min(uxm,0).^2 + max(uxp,0).^2 + min(uym,0).^2 + max(uyp,0).^2).^0.5;
    % K*du     
    k = (ux.^2.*uyy + uy.^2.*uxx - 2*ux.*uy.*uxy)./max(du.^2 , 1e-10); 
    % u
    u = u + dt*(g.*k + alpha*max(g,0).*dpu + alpha*min(g,0).*dmu + ...
       max(gx,0).*uxm + min(gx,0).*uxp + max(gy,0).*uym + min(gy,0).*uyp);

    if mod(i, Reinit_step) == 0
      u = Reinitial2D(u, Reinit_step);
    end 
    if mod(i, 10) == 0
      pause(0.001);
      imshow(I,[]);
      hold on;
      M = max(max(u));
      contour(u, [-tol*M tol*M], 'r');
      hold off;
    end

end 
u1 = u;
M = max(max(u));
figure
fu = abs(u > 1e-16 * M);
imshow(fu,[])

figure 
imshow(g,[]);
end 
