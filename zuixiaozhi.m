function y=zuixiaozhi(f,g,h,m)

ymin=h;

for t=0.5/10000:0.5/10000:0.5
    if ymin>f*t^3+g*t^2+h*t+m
        ymin=f*t^3+g*t^2+h*t+m;
    end
end
y=ymin;
end



