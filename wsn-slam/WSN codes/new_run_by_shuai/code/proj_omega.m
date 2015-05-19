function y=proj_omega(lower,upper,x);
y=x.*(lower<=x).*(upper>=x)+lower.*(lower>x)+upper.*(upper<x);
if (lower==-inf)&(upper==inf)
    y=x;
end