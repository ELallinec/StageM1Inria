function Ha_tab = parareal_Ha(T,y0, dtf,dtg,N)
  delta_t = T/N;
  kmax = 6;
  Ha_tab = zeros(kmax+2, N+1);
  lambda_tab_old = zeros (2, N+1);
  lambda_tab = zeros (2, N+1);
  lambda_f = zeros (2, N+1);
  lambda_g = zeros (2, N+1);
  tab_f = zeros(2,N+1);
  lambda_tab_old(:,1)= y0;
  lambda_tab(:,1)= y0;
  lambda_f(:,1)= y0;
  lambda_g(:,1)= y0;
  tab_f(:,1)=y0;
  for n=2:N+1
    tab_f(:,n) = Stormer((n-2)*delta_t,(n-1)*delta_t,tab_f(:,n-1),dtf);
  end
  Ha_tab(1,:) = 0.5*(tab_f(1,:).^2 + tab_f(2,:).^2);
  for n=2:N+1
    lambda_tab_old(:,n) = Stormer((n-2)*delta_t,(n-1)*delta_t,lambda_tab_old(:,n-1),dtg);
  end
  Ha_tab(2,:) = 0.5*(lambda_tab_old(1,:).^2 + lambda_tab_old(2,:).^2);
  k=1;
  while(k<=kmax)  
    for n=2:N+1
      lambda_f(:,n) = Stormer((n-2)*delta_t,(n-1)*delta_t,lambda_tab_old(:,n-1),dtf);
      lambda_g(:,n)= Stormer((n-2)*delta_t,(n-1)*delta_t,lambda_tab_old(:,n-1),dtg);
    end
    for n =2:N+1 
      lambda_tab(:,n) = lambda_f(:,n) - lambda_g(:,n) + Stormer((n-2)*delta_t,(n-1)*delta_t,lambda_tab(:,n-1),dtg);
      
    end
    Ha_tab(k+2,:) = 0.5*(lambda_tab(1,:).^2 + lambda_tab(2,:).^2);
    lambda_tab_old = lambda_tab;
    disp(k);
    k +=1;
  end
end