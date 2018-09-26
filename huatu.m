syms t;

apha_final=@(t) y1(33,7)*t.^3+y2(33,7)*t.^2+y3(33,7)*t+y4(33,7);
sita_final=@(t) y5(33,7)*t.^3+y6(33,7)*t.^2+y7(33,7)*t+y8(33,7);

t=0:0.5/100:0.5;

plot(t,-0.9*cos(apha_final(t))+0.9*cos(apha_final(t)+sita_final(t)));  %TOE的x坐标
title('TOE x POSITION');

plot(t,0.9*sin(sita_final(t))-0.9*sin(sita_final(t)+apha_final(t)));  %TOE的y坐标
title('TOE y POSITION');

plot(t,0.9*sin(apha_final)*diff(apha_final)-0.5*0.9*sin(sita_final+apha_final)*diff(sita_final+apha_final));  %TOE的x速度
title('TOE x VELOCITY');

plot(t,0.9*cos(apha_final)*diff(apha_final)-0.5*0.9*cos(sita_final+apha_final)*diff(sita_final+apha_final));  %TOE的y速度
title('TOE y VELOCITY');

plot(-0.9*cos(apha_final(t))+0.9*cos(apha_final(t)+sita_final(t)),0.9*sin(sita_final(t))-0.9*sin(sita_final(t)+apha_final(t)));  %TOE的轨迹
title('TOE TRAJECTORY');

plot(t,-0.9*cos(apha_final));  %HIP的x坐标
title('HIP x POSITION');

plot(t,0.9*sin(apha_final));  %HIP的y坐标
title('HIP y POSITION');

plot(t,diff(-0.9*cos(apha_final)));  %HIP的x速度
title('HIP x VELOCITY');

plot(t,diff(0.9*sin(apha_final)));  %HIP的y速度
title('HIP y VELOCITY');

plot(-0.9*cos(apha_final),0.9*sin(apha_final));  %HIP的轨迹
title('HIP TRAJECTORY');

plot(t,M1_final);  %TOE的M力矩
title('TOE Moment');

plot(t,M2_final);  %HIP的M力矩
title('HIP Moment');