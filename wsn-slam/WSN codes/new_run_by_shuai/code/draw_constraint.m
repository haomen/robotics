function draw_constraint(x,d,eta_plus,eta_minus,theta,delta,fig_handle)
d0=d+eta_plus;
d1=d-eta_minus;
if d1<0
    d1=0;
end

angles=[-delta:0.01:delta,delta];
theta=(theta+angles);
x01=x(1)+d0*cos(theta);
x02=x(2)+d0*sin(theta);
x11=x(1)+d1*cos(theta);
x12=x(2)+d1*sin(theta);
figure(fig_handle),hold on,
plot(x01,x02);
plot(x11,x12);
plot([x01(1) x11(1)],[x02(1),x12(1)]);
plot([x01(end) x11(end)],[x02(end),x12(end)]);