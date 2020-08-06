z = [-100:0.001:100];


bar = zeros(1,length(z));
for i = 1:length(z)
    bar(i) = logbarrierfunc(0.5,z(i));
end

plot(z,bar);