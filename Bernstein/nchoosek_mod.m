function out = nchoosek_mod(n,k)

out = 1;
for i = 1:k
    out = out*(n-(k-i));
    out = out/i;
end

end