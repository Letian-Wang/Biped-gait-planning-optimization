load data;
k=[0];
m=[0];
n=[0];
for i=1:201
    for j=1:201
          
       
        %apha_=@(t) y1(i,j)*t^3+y2(i,j)*t^2+y3(i,j)*t+y4(i,j);
        %sita_=@(t) y5(i,j)*t^3+y6(i,j)*t^2+y7(i,j)*t+y8(i,j); %为了满足function handle
        
        
        
        %[x,fval]=fminbnd(@(t) 0.9*sin(apha_(t))-0.9*sin(apha_(t)+sita_(t)),0,0.5); %最低点
        %  if fval>0
     
        if zuixiaozhi(y1(i,j),y2(i,j),y3(i,j),y4(i,j),y5(i,j),y6(i,j),y7(i,j),y8(i,j))>0
            m=[m i];
            n=[n j];
             k=[k zuixiaozhi(y1(i,j),y2(i,j),y3(i,j),y4(i,j),y5(i,j),y6(i,j),y7(i,j),y8(i,j))];
        end
            
       
    end
end
save data2.mat；
%apha_=@(t) -0.94*t^3+26734*t^2-13364*t+0.5553;
 %   sita_=@(t) -0.68*t^3+6235.3*t^2-3121.5*t+1.1262;

  %  t=0:0.5/100:0.5;
%plot(t,0.9*sin(-0.94.*t.^3+26734.*t.^2-13364.*t+0.5553)-0.9.*sin(-0.94.*t.^3+26734.*t.^2-13364.*t+0.5553+-0.68.*t.^3+6235.3.*t.^2-3121.5.*t+1.1262))
%