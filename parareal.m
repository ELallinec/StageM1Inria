function lambda_tab = parareal(T,y0, dtf,dtg,N)
  delta_t = T/N;
  kmax = 10;
  lambda_tab_old = zeros (2, N+1);
  lambda_tab = zeros (2, N+1);
  lambda_f = zeros (2, N+1);
  lambda_g = zeros (2, N+1)
  lambda_tab_old(:,1)= y0;
  lambda_tab(:,1)= y0;
  for n=2:N+1
    lambda_tab_old(:,n) = Stormer((n-2)*delta_t,(n-1)*delta_t,lambda_tab_old(:,n-1),dtg);
  endfor
  k=1;
  while(k<=kmax)
    lambda_f(:,1)= y0;
    lambda_g(:,1)= y0;
    for n=2:N+1
      lambda_f(:,n) = Stormer((n-2)*delta_t,(n-1)*delta_t,lambda_tab_old(:,n-1),dtf);
      lambda_g(:,n)= Stormer((n-2)*delta_t,(n-1)*delta_t,lambda_tab_old(:,n-1),dtg);
      endfor
    for n =2:N+1 
      lambda_tab(:,n) = lambda_f(:,n) - lambda_g(:,n) + Stormer((n-2)*delta_t,(n-1)*delta_t,lambda_tab(:,n-1),dtg);
    endfor
    lambda_tab_old = lambda_tab;
    disp(k);
    k +=1;
  endwhile
endfunction