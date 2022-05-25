T = 10;
N = 10*T;
delta_T = T/N;
dtg = delta_T;
dtf= dtg/100;
y0 = [0; 1];

%y =parareal(T,y0,dtf,dtg,N);
%plot(y(2,:),y(1,:));
Ha = parareal_Ha(T,y0,dtf,dtg,N);
tab_t = linspace(0,T,N+1);
H0 = 0.5*ones(1,N+1);
for k=1:rows(Ha)
  loglog(tab_t(2:end),norm(Ha(k,2:end)-H0(2:end)),'-');
  hold on;
endfor

