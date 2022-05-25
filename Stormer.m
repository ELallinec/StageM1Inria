function y = Stormer(t0,T,y0,dt)
	m = (T-t0)/dt;
	p = y0(1);
	q = y0(2);
  y = y0;
	for i=1:m+1
	  q_temp = q + (dt / 2) * p;
    p = p - dt * q_temp;
    q = q_temp +(dt / 2) * p;
    
	endfor
  y = [p;q];
endfunction