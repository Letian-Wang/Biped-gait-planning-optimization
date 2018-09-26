function y=zuixiaozhi2(e,f,g,h)

ymin=0.125*e+0.25*f+0.5*g+h;

for t=0.15:0.35/10000:0.5
    if ymin>e*t^3+f*t^2+g*t+h
        ymin=e*t^3+f*t^2+g*t+h;
    end
end
y=ymin;
end

