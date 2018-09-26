syms t;

apha_final=@(t) y1(33,7)*t.^3+y2(33,7)*t.^2+y3(33,7)*t+y4(33,7);
sita_final=@(t) y5(33,7)*t.^3+y6(33,7)*t.^2+y7(33,7)*t+y8(33,7);

t=0:0.5/100:0.5;

plot(t,-0.9*cos(apha_final(t))+0.9*cos(apha_final(t)+sita_final(t)));  %TOE��x����
title('TOE x POSITION');

plot(t,0.9*sin(sita_final(t))-0.9*sin(sita_final(t)+apha_final(t)));  %TOE��y����
title('TOE y POSITION');

plot(t,0.9*sin(apha_final)*diff(apha_final)-0.5*0.9*sin(sita_final+apha_final)*diff(sita_final+apha_final));  %TOE��x�ٶ�
title('TOE x VELOCITY');

plot(t,0.9*cos(apha_final)*diff(apha_final)-0.5*0.9*cos(sita_final+apha_final)*diff(sita_final+apha_final));  %TOE��y�ٶ�
title('TOE y VELOCITY');

plot(-0.9*cos(apha_final(t))+0.9*cos(apha_final(t)+sita_final(t)),0.9*sin(sita_final(t))-0.9*sin(sita_final(t)+apha_final(t)));  %TOE�Ĺ켣
title('TOE TRAJECTORY');

plot(t,-0.9*cos(apha_final));  %HIP��x����
title('HIP x POSITION');

plot(t,0.9*sin(apha_final));  %HIP��y����
title('HIP y POSITION');

plot(t,diff(-0.9*cos(apha_final)));  %HIP��x�ٶ�
title('HIP x VELOCITY');

plot(t,diff(0.9*sin(apha_final)));  %HIP��y�ٶ�
title('HIP y VELOCITY');

plot(-0.9*cos(apha_final),0.9*sin(apha_final));  %HIP�Ĺ켣
title('HIP TRAJECTORY');

plot(t,M1_final);  %TOE��M����
title('TOE Moment');

plot(t,M2_final);  %HIP��M����
title('HIP Moment');