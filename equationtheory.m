a = @(t) t;
w = 10*rand(1);


eq

b(t)=a(t)+cos(wt-k[d-b(t)])




function [c,ceq] =  nonlcon(coefs,a,w,k,d)
    bigeq = @(t) b(coefs,t) - a(t) - cos(w*t-k*(d-b(coefs,t)));
    c = integral(bigeq,tmin,tmax);
    ceq = [];
end