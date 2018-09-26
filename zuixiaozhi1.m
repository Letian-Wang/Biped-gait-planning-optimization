function y=zuixiaozhi1(e,f,g,h)

ymin=h;

for t=0.15/10000:0.15/10000:0.15
    if ymin>e*t^3+f*t^2+g*t+h
        ymin=e*t^3+f*t^2+g*t+h;
    end
end
y=ymin;
end

