function tab_y = tab_Stormer(t0,T,y0,dt)
	m = (T-t0)/dt;
	p = y0(1);
	q = y0(2);
  tab_y = zeros(2,m+1);
	for i=1:m+1
    tab_y(:,i) = [p;q]
	  p_temp = p - (dt / 2) * q;
    q = q + dt * p_temp;
    p = p_temp -(dt / 2) * q;
	endfor
 
endfunction