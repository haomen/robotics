function y=map(theta)
theta1=mod(theta,2*pi);
if theta1>=pi
    theta1=theta-2*pi;
end
y=(theta1>=pi).*(theta1-2*pi)+(theta1<pi).*(theta1);
