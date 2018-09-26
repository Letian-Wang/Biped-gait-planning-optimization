syms t;
j=33;
i=7;
m=7;
n=33;
 apha=y1(i,j)*t^3+y2(i,j)*t^2+y3(i,j)*t+y4(i,j);
        sita=y5(i,j)*t^3+y6(i,j)*t^2+y7(i,j)*t+y8(i,j); %变化中的
        
        M11=-1/2*0.9*1*(-1/2*0.9*(diff(apha,2)+diff(sita,2))*(cos(apha+sita))^2+(0.9*(diff(sita))^2*sin(apha)+2*0.9*diff(apha)*diff(sita)*sin(apha)+2*0.9*cos(apha)*diff(apha,2)+0.9*diff(sita,2)*cos(apha)+1)*cos(apha+sita)-1/2*0.9*(diff(apha,2)+diff(sita,2))*(sin(apha+sita))^2-0.9*(-2*sin(apha)*diff(apha,2)-sin(apha)*diff(sita,2)+diff(sita)*cos(apha)*(diff(sita)+2*diff(apha)))*sin(apha+sita)-5/2*0.9*((cos(apha))^2+(sin(apha))^2+1/15)*diff(apha,2)-cos(apha)*(9.8*2));
        M22=1/4*0.9*(0.9*(diff(apha,2)+diff(sita,2))*(cos(apha+sita))^2+(2*0.9*sin(apha)*(diff(apha))^2-2*0.9*cos(apha)*diff(apha,2)-2)*cos(sita+apha)+0.9*((diff(apha,2)+diff(sita,2))*(sin(apha+sita))^2+(-2*cos(apha)*(diff(apha,2))^2-2*sin(apha)*diff(apha,2))*sin(apha+sita)+1/3*diff(sita,2)));
        
  
        
        
        Emin=int(M11,t,0,0.5)+int(M22,t,0,0.5);  %消耗能量
for j=84:89
        apha=y1(1,j)*t^3+y2(1,j)*t^2+y3(1,j)*t+y4(1,j);
        sita=y5(1,j)*t^3+y6(1,j)*t^2+y7(1,j)*t+y8(1,j); %变化中的
        
        M11=-1/2*0.9*1*(-1/2*0.9*(diff(apha,2)+diff(sita,2))*(cos(apha+sita))^2+(0.9*(diff(sita))^2*sin(apha)+2*0.9*diff(apha)*diff(sita)*sin(apha)+2*0.9*cos(apha)*diff(apha,2)+0.9*diff(sita,2)*cos(apha)+1)*cos(apha+sita)-1/2*0.9*(diff(apha,2)+diff(sita,2))*(sin(apha+sita))^2-0.9*(-2*sin(apha)*diff(apha,2)-sin(apha)*diff(sita,2)+diff(sita)*cos(apha)*(diff(sita)+2*diff(apha)))*sin(apha+sita)-5/2*0.9*((cos(apha))^2+(sin(apha))^2+1/15)*diff(apha,2)-cos(apha)*(9.8*2));
        M22=1/4*0.9*(0.9*(diff(apha,2)+diff(sita,2))*(cos(apha+sita))^2+(2*0.9*sin(apha)*(diff(apha))^2-2*0.9*cos(apha)*diff(apha,2)-2)*cos(sita+apha)+0.9*((diff(apha,2)+diff(sita,2))*(sin(apha+sita))^2+(-2*cos(apha)*(diff(apha,2))^2-2*sin(apha)*diff(apha,2))*sin(apha+sita)+1/3*diff(sita,2)));
        

       
        E=int(M11,t,0,0.5)+int(M22,t,0,0.5);  %消耗能量
        if E<Emin
            Emin=E; %最低能量
            m=1;
            n=j;
        end
end
            
            