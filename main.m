tic
%%%%%计算不同取值条件下的变量
n=-1:0.01:1;
m=-1:0.01:1;
y1=zeros(201,201);
y2=zeros(201,201);
y3=zeros(201,201);
y4=zeros(201,201);
y5=zeros(201,201);
y6=zeros(201,201);
y7=zeros(201,201);
y8=zeros(201,201);
y9=zeros(201,201);
y10=zeros(201,201);
y11=zeros(201,201);
y12=zeros(201,201);
y13=zeros(201,201);
y14=zeros(201,201);
y15=zeros(201,201);
y16=zeros(201,201);

F0=[0;0;0;0;0;0;0;0;0;0;0;0];
options=optimset('Display','off');
for i=1:201
for j=1:201
i
j
F0=fsolve(@(x)qiujie_shijie(x,n(i),m(j)),F0,options);
y1(i,j)=n(i);
y2(i,j)=F0(1);
y3(i,j)=F0(2);
y4(i,j)=F0(3);
y5(i,j)=m(j);
y6(i,j)=F0(4);
y7(i,j)=F0(5);
y8(i,j)=F0(6);
y9(i,j)=n(i);
y10(i,j)=F0(7);
y11(i,j)=F0(8);
y12(i,j)=F0(9);
y13(i,j)=m(j);
y14(i,j)=F0(10);
y15(i,j)=F0(11);
y16(i,j)=F0(12);
end
end
toc

%tihuan=y5;
%y5=y6;
%y6=tihuan;

%%%%%筛选出不会触地的变量组合
gaoduzuidi=zeros(201,201);
ayouxiao=[0];
eyouxiao=[0];

for i=1:201
    for j=1:201    
        i
        j
        
        gaoduzuidi1(i,j)=zuixiaozhi1(y5(i,j),y6(i,j),y7(i,j),y8(i,j));
        gaoduzuidi2(i,j)=zuixiaozhi2(y13(i,j),y14(i,j),y15(i,j),y16(i,j));
        
        
        %gaoduzuidi1(i,j)=zuixiaozhi(y5(i,j),y6(i,j),y7(i,j),y8(i,j));
        %gaoduzuidi2(i,j)=zuixiaozhi(y13(i,j),y14(i,j),y15(i,j),y16(i,j));
        if gaoduzuidi1(i,j)>=-0.1&&gaoduzuidi2(i,j)>=-0.1
            
           ayouxiao=[ayouxiao i];
           eyouxiao=[eyouxiao j];
        end
    end
end
save data2.mat




%%%能量最低



tic
syms t;
i=2

        %x_=y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i));
        %y_=y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i));
        
        %奇异前的计算
        shijian_qian=0.15;%避免奇异性
       
        x_=y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i));
        y_=y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i));
        
        sita_qian=acos((2*0.8^2-(x_)^2-(y_)^2)/(2*0.8^2));
        apha_qian=acos(((x_)^2+(y_)^2)/(2*0.8^2))+atan((y_)/(-(x_)));
        
        
        M11=-1/2*0.9*1*(-1/2*0.9*(diff(apha_qian,2)+diff(sita_qian,2))*(cos(apha_qian+sita_qian))^2+(0.9*(diff(sita_qian))^2*sin(apha_qian)+2*0.9*diff(apha_qian)*diff(sita_qian)*sin(apha_qian)+2*0.9*cos(apha_qian)*diff(apha_qian,2)+0.9*diff(sita_qian,2)*cos(apha_qian)+1)*cos(apha_qian+sita_qian)-1/2*0.9*(diff(apha_qian,2)+diff(sita_qian,2))*(sin(apha_qian+sita_qian))^2-0.9*(-2*sin(apha_qian)*diff(apha_qian,2)-sin(apha_qian)*diff(sita_qian,2)+diff(sita_qian)*cos(apha_qian)*(diff(sita_qian)+2*diff(apha_qian)))*sin(apha_qian+sita_qian)-5/2*0.9*((cos(apha_qian))^2+(sin(apha_qian))^2+1/15)*diff(apha_qian,2)-cos(apha_qian)*(9.8*2));
        M12=1/4*0.9*(0.9*(diff(apha_qian,2)+diff(sita_qian,2))*(cos(apha_qian+sita_qian))^2+(2*0.9*sin(apha_qian)*(diff(apha_qian))^2-2*0.9*cos(apha_qian)*diff(apha_qian,2)-2)*cos(sita_qian+apha_qian)+0.9*((diff(apha_qian,2)+diff(sita_qian,2))*(sin(apha_qian+sita_qian))^2+(-2*cos(apha_qian)*(diff(apha_qian,2))^2-2*sin(apha_qian)*diff(apha_qian,2))*sin(apha_qian+sita_qian)+1/3*diff(sita_qian,2)));
        
        sita_nengliang_1=sita_qian*M12;
        apha_nengliang_1=apha_qian*M11;

        E1=int(apha_nengliang_1,t,0,shijian_qian)+int(sita_nengliang_1,t,0,shijian_qian);  %消耗能量
 
        
        
        
        %奇异后的计算
        shijian_hou=0.15;
        
        x_=y9(ayouxiao(i),eyouxiao(i))*t^3+y10(ayouxiao(i),eyouxiao(i))*t^2+y11(ayouxiao(i),eyouxiao(i))*t+y12(ayouxiao(i),eyouxiao(i));
        y_=y13(ayouxiao(i),eyouxiao(i))*t^3+y14(ayouxiao(i),eyouxiao(i))*t^2+y15(ayouxiao(i),eyouxiao(i))*t+y16(ayouxiao(i),eyouxiao(i));
        
        
        
        sita_hou=-acos((2*0.8^2-(x_)^2-(y_)^2)/(2*0.8^2));
        apha_hou=pi-(acos(((x_)^2+(y_)^2)/(2*0.8^2))+atan((y_)/(x_)));
        
        M21=-1/2*0.9*1*(-1/2*0.9*(diff(apha_hou,2)+diff(sita_hou,2))*(cos(apha_hou+sita_hou))^2+(0.9*(diff(sita_hou))^2*sin(apha_hou)+2*0.9*diff(apha_hou)*diff(sita_hou)*sin(apha_hou)+2*0.9*cos(apha_hou)*diff(apha_hou,2)+0.9*diff(sita_hou,2)*cos(apha_hou)+1)*cos(apha_hou+sita_hou)-1/2*0.9*(diff(apha_hou,2)+diff(sita_hou,2))*(sin(apha_hou+sita_hou))^2-0.9*(-2*sin(apha_hou)*diff(apha_hou,2)-sin(apha_hou)*diff(sita_hou,2)+diff(sita_hou)*cos(apha_hou)*(diff(sita_hou)+2*diff(apha_hou)))*sin(apha_hou+sita_hou)-5/2*0.9*((cos(apha_hou))^2+(sin(apha_hou))^2+1/15)*diff(apha_hou,2)-cos(apha_hou)*(9.8*2));
        M22=1/4*0.9*(0.9*(diff(apha_hou,2)+diff(sita_hou,2))*(cos(apha_hou+sita_hou))^2+(2*0.9*sin(apha_hou)*(diff(apha_hou))^2-2*0.9*cos(apha_hou)*diff(apha_hou,2)-2)*cos(sita_hou+apha_hou)+0.9*((diff(apha_hou,2)+diff(sita_hou,2))*(sin(apha_hou+sita_hou))^2+(-2*cos(apha_hou)*(diff(apha_hou,2))^2-2*sin(apha_hou)*diff(apha_hou,2))*sin(apha_hou+sita_hou)+1/3*diff(sita_hou,2)));
        
        sita_nengliang_2=sita_qian*M22;
        apha_nengliang_2=apha_qian*M21;
        
        E2=int(apha_nengliang_2,t,shijian_hou,0.5)+int(sita_nengliang_2,t,shijian_hou,0.5);  %消耗能量
        
     
        E=E1+E2;
 
toc



tic
for i=2:15220
    cunzai=0;
        i
        %求解发生奇异性的时间
        jie=solve(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)),t);
        k=length(jie);
        for j=1:k
            tf=isreal(jie(j));
            if tf==1&&jie(j)>0&&jie(j)<0.5
                shijian=jie(j);
                cunzai=1;
            end
        end
        if cunzai==1
             x_=y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i));
        y_=y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i));
        
        
        
        
        %奇异前的计算
        shijian_qian=shijian-0.001;%避免奇异性
        sita_qian=acos((2*0.8^2-(x_)^2-(y_)^2)/(2*0.8^2));
        apha_qian=acos(((x_)^2+(y_)^2)/(2*0.8^2))+atan((y_)/(-(x_)));
        
        M11=-1/2*0.9*1*(-1/2*0.9*(diff(apha_qian,2)+diff(sita_qian,2))*(cos(apha_qian+sita_qian))^2+(0.9*(diff(sita_qian))^2*sin(apha_qian)+2*0.9*diff(apha_qian)*diff(sita_qian)*sin(apha_qian)+2*0.9*cos(apha_qian)*diff(apha_qian,2)+0.9*diff(sita_qian,2)*cos(apha_qian)+1)*cos(apha_qian+sita_qian)-1/2*0.9*(diff(apha_qian,2)+diff(sita_qian,2))*(sin(apha_qian+sita_qian))^2-0.9*(-2*sin(apha_qian)*diff(apha_qian,2)-sin(apha_qian)*diff(sita_qian,2)+diff(sita_qian)*cos(apha_qian)*(diff(sita_qian)+2*diff(apha_qian)))*sin(apha_qian+sita_qian)-5/2*0.9*((cos(apha_qian))^2+(sin(apha_qian))^2+1/15)*diff(apha_qian,2)-cos(apha_qian)*(9.8*2));
        M12=1/4*0.9*(0.9*(diff(apha_qian,2)+diff(sita_qian,2))*(cos(apha_qian+sita_qian))^2+(2*0.9*sin(apha_qian)*(diff(apha_qian))^2-2*0.9*cos(apha_qian)*diff(apha_qian,2)-2)*cos(sita_qian+apha_qian)+0.9*((diff(apha_qian,2)+diff(sita_qian,2))*(sin(apha_qian+sita_qian))^2+(-2*cos(apha_qian)*(diff(apha_qian,2))^2-2*sin(apha_qian)*diff(apha_qian,2))*sin(apha_qian+sita_qian)+1/3*diff(sita_qian,2)));
        
        sita_nengliang_1=sita_qian*M12;
        apha_nengliang_1=apha_qian*M11;

        E1=int(apha_nengliang_1,t,0,shijian_qian)+int(sita_nengliang_1,t,0,shijian_qian);  %消耗能量
        %sita_qian=acos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2));
        %apha_qian=acos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))));

        
        
        %y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i))
         %sita=arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))
        %apha=arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))));

        %sita=arccos((2*0.8^2-(a*t^3+b*t^2+c*t+d)^2-(e*t^3+f*t^2+g*t+h)^2)/(2*0.8^2))
        %apha=arcos(((a*t^3+b*t^2+c*t+d)^2+(e*t^3+f*t^2+g*t+h)^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(-(a*t^3+b*t^2+c*t+d)));

         
       %%% M11=-1/2*0.9*1*(-1/2*0.9*(diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))),2)+diff(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2))*(cos(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))))+arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))))^2+(0.9*(diff(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))))^2*sin(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))+2*0.9*diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))*diff(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)))*sin(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))+2*0.9*cos(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))*diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))),2)+0.9*diff(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2)*cos(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))+1)*cos(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))))+arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)))-1/2*0.9*(diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))),2)+diff(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2))*(sin(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))))+arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))))^2-0.9*(-2*sin(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))*diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))),2)-sin(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))*diff(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2)+diff(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)))*cos(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))*(diff(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)))+2*diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))))*sin(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))))+arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)))-5/2*0.9*((cos(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))))))^2+(sin(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))))))^2+1/15)*diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))),2)-cos(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))*(9.8*2));
        %%%M12=1/4*0.9*(0.9*(diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))),2)+diff(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2))*(cos(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))))+arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))))^2+(2*0.9*sin(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))*(diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))))))^2-2*0.9*cos(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))*diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))),2)-2)*cos(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))+0.9*((diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))),2)+diff(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2))*(sin(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))))+arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))))^2+(-2*cos(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))*(diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))),2))^2-2*sin(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))))*diff(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))),2))*sin(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))))+arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)))+1/3*diff(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2)));
        %%%apha_nengliang_1=M11*(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2+( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan(( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))/(-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i)))));
        %%%sita_nengliang_1=M12*(arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i))^2-( y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)));
        
        
        
        
        
        %奇异后的计算
        shijian_hou=shijian+0.001;%避免奇异性
        sita_hou=-acos((2*0.8^2-(x_)^2-(y_)^2)/(2*0.8^2));
        apha_hou=pi-(acos(((x_)^2+(y_)^2)/(2*0.8^2))+atan((y_)/(x_)));
        
        M21=-1/2*0.9*1*(-1/2*0.9*(diff(apha_hou,2)+diff(sita_hou,2))*(cos(apha_hou+sita_hou))^2+(0.9*(diff(sita_hou))^2*sin(apha_hou)+2*0.9*diff(apha_hou)*diff(sita_hou)*sin(apha_hou)+2*0.9*cos(apha_hou)*diff(apha_hou,2)+0.9*diff(sita_hou,2)*cos(apha_hou)+1)*cos(apha_hou+sita_hou)-1/2*0.9*(diff(apha_hou,2)+diff(sita_hou,2))*(sin(apha_hou+sita_hou))^2-0.9*(-2*sin(apha_hou)*diff(apha_hou,2)-sin(apha_hou)*diff(sita_hou,2)+diff(sita_hou)*cos(apha_hou)*(diff(sita_hou)+2*diff(apha_hou)))*sin(apha_hou+sita_hou)-5/2*0.9*((cos(apha_hou))^2+(sin(apha_hou))^2+1/15)*diff(apha_hou,2)-cos(apha_hou)*(9.8*2));
        M22=1/4*0.9*(0.9*(diff(apha_hou,2)+diff(sita_hou,2))*(cos(apha_hou+sita_hou))^2+(2*0.9*sin(apha_hou)*(diff(apha_hou))^2-2*0.9*cos(apha_hou)*diff(apha_hou,2)-2)*cos(sita_hou+apha_hou)+0.9*((diff(apha_hou,2)+diff(sita_hou,2))*(sin(apha_hou+sita_hou))^2+(-2*cos(apha_hou)*(diff(apha_hou,2))^2-2*sin(apha_hou)*diff(apha_hou,2))*sin(apha_hou+sita_hou)+1/3*diff(sita_hou,2)));
        
        sita_nengliang_2=sita_qian*M22;
        apha_nengliang_2=apha_qian*M21;
        
        E2=int(apha_nengliang_2,t,shijian_hou,0.5)+int(sita_nengliang_2,t,shijian_hou,0.5);  %消耗能量
        
        %sita=-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2));
        %apha=pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))));
        
        
        
        %sita=-arccos((2*0.8^2-(a*t^3+b*t^2+c*t+d)^2-(e*t^3+f*t^2+g*t+h)^2)/(2*0.8^2));
        %apha=pi-(arcos(((a*t^3+b*t^2+c*t+d)^2+(e*t^3+f*t^2+g*t+h)^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(a*t^3+b*t^2+c*t+d)));
        
        %M21=-1/2*0.9*1*(-1/2*0.9*(diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))),2)+diff(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2))*(cos(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))))+-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))))^2+(0.9*(diff(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))))^2*sin(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))+2*0.9*diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))*diff(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)))*sin(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))+2*0.9*cos(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))*diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))),2)+0.9*diff(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2)*cos(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))+1)*cos(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))))+-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)))-1/2*0.9*(diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))),2)+diff(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2))*(sin(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))))+-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))))^2-0.9*(-2*sin(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))*diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))),2)-sin(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))*diff(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2)+diff(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)))*cos(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))*(diff(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)))+2*diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))))*sin(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))))+-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)))-5/2*0.9*((cos(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))))))^2+(sin(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))))))^2+1/15)*diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))),2)-cos(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))*(9.8*2));
        %M22=1/4*0.9*(0.9*(diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))),2)+diff(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2))*(cos(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))))+-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))))^2+(2*0.9*sin(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))*(diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))))))^2-2*0.9*cos(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))*diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))),2)-2)*cos(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))+0.9*((diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))),2)+diff(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2))*(sin(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))))+-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))))^2+(-2*cos(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))*(diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))),2))^2-2*sin(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))))*diff(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))),2))*sin(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))))+-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)))+1/3*diff(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)),2)));
        %apha_nengliang_2=M21*(pi-(arcos(((y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2+(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2))+arctan((e*t^3+f*t^2+c*t+d)/(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i))))));
        %sita_nengliang_2=M22*(-arccos((2*0.8^2-(y1(ayouxiao(i),eyouxiao(i))*t^3+y2(ayouxiao(i),eyouxiao(i))*t^2+y3(ayouxiao(i),eyouxiao(i))*t+y4(ayouxiao(i),eyouxiao(i)))^2-(y5(ayouxiao(i),eyouxiao(i))*t^3+y6(ayouxiao(i),eyouxiao(i))*t^2+y7(ayouxiao(i),eyouxiao(i))*t+y8(ayouxiao(i),eyouxiao(i)))^2)/(2*0.8^2)));
        
      
        %E2=int(apha_nengliang_2,t,shijian,0.5)+int(sita_nengliang_2,t,shijian,0.5);  %消耗能量
        
        
        %%总的能量
        E=E1+E2;
        if E<Emin&&E>0
            Emin=E;
            i_zuixiao=i;
        end
        end
            
       
end
toc
save data3.mat            
       
%%画图
i=1;
j=1;

%solve(x_hip^2+y_hip^2==0.8^2,(x_hip-(y1(i,j)*t^3+y2(i,j)*t^2+y3(i,j)*t+y4(i,j)))^2+(y_hip-(y5(i,j)*t^3+y6(i,j)*t^2+y7(i,j)*t+y8(i,j)))^==0.8^2,


t=0:0.5/1000:0.5;

for s=1:length(t)
    if t(s)<0.15
        x__=y1(i,j)*t.^3+y2(i,j)*t.^2+y3(i,j)*t+y4(i,j);
        y__=y5(i,j)*t.^3+y6(i,j)*t.^2+y7(i,j)*t+y8(i,j);
    else
        x__=y9(i,j)*t.^3+y10(i,j)*t.^2+y11(i,j)*t+y12(i,j);
        y__=y13(i,j)*t.^3+y14(i,j)*t.^2+y15(i,j)*t+y16(i,j);
    end
end
figure
subplot(1,2,1)
plot(t,x__)
title('TOE x POSITION');
subplot(1,2,2)
plot(t,y__)
title('TOE y POSITION')
    


t=0:0.15/1000:0.15;
 
figure
subplot(1,2,1)
%x_=y1(i,j)*t^3+y2(i,j)*t^2+y3(i,j)*t+y4(i,j);
%y_=y5(i,j)*t^3+y6(i,j)*t^2+y7(i,j)*t+y8(i,j); 
plot(t,y1(i,j)*t.^3+y2(i,j)*t.^2+y3(i,j)*t+y4(i,j));  %TOE的x坐标
title('TOE x POSITION');
grid on
subplot(1,2,2)
plot(t,y5(i,j)*t.^3+y6(i,j)*t.^2+y7(i,j)*t+y8(i,j));  %TOE的y坐标
title('TOE y POSITION');
grid on

figure
subplot(1,2,1)
plot(t,3*y1(i,j)*t.^2+2*y2(i,j)*t+y3(i,j));  %TOE的x速度
title('TOE x VELOCITY');
subplot(1,2,2)
plot(t,3*y5(i,j)*t.^2+2*y6(i,j)*t+y7(i,j));  %TOE的y速度
title('TOE y VELOCITY');

figure
plot(y1(i,j)*t.^3+y2(i,j)*t.^2+y3(i,j)*t+y4(i,j),y5(i,j)*t.^3+y6(i,j)*t.^2+y7(i,j)*t+y8(i,j));  %TOE的轨迹
title('TOE TRAJECTORY');

t=0.15:0.35/100:0.35;
 








figure
subplot(1,2,1)
%x_=y1(i,j)*t^3+y2(i,j)*t^2+y3(i,j)*t+y4(i,j);
%y_=y5(i,j)*t^3+y6(i,j)*t^2+y7(i,j)*t+y8(i,j); 
plot(t,y9(i,j)*t.^3+y10(i,j)*t.^2+y11(i,j)*t+y12(i,j));  %TOE的x坐标
title('TOE x POSITION');
grid on
subplot(1,2,2)
plot(t,y13(i,j)*t.^3+y14(i,j)*t.^2+y15(i,j)*t+y16(i,j));  %TOE的y坐标
title('TOE y POSITION');
grid on

figure
subplot(1,2,1)
plot(t,3*y9(i,j)*t.^2+2*y10(i,j)*t+y11(i,j));  %TOE的x速度
title('TOE x VELOCITY');
subplot(1,2,2)
plot(t,3*y12(i,j)*t.^2+2*y13(i,j)*t+y14(i,j));  %TOE的y速度
title('TOE y VELOCITY');

figure
plot(y9(i,j)*t.^3+y10(i,j)*t.^2+y11(i,j)*t+y12(i,j),y13(i,j)*t.^3+y14(i,j)*t.^2+y15(i,j)*t+y16(i,j));  %TOE的轨迹
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
