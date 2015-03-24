% localization in wsn via a distributed dual neural network
clear all
close all
m=9;%number of anchors
n=150;%number of blind nodes
n=m+n;%number of total sensors
Error_level=0.00;
%geometry graph
%cretical distance;
Cretical_dis=0.15;
Coor_x=rand(n,1);
Coor_y=rand(n,1);

%Coor_x(1:m)=round(Coor_x(1:m));
%Coor_y(1:m)=round(Coor_y(1:m));
mmm=sqrt(m);
x111=0:1/(mmm-1):1;
y111=x111;
[xy111,yx111]=meshgrid(x111,y111);


Coor_x(1:m)=xy111(1:end);
Coor_y(1:m)=yx111(1:end);

clear i;
Coor=Coor_x+i*Coor_y;
Dis_matrix=Coor*ones(1,n)-ones(n,1)*conj(Coor');
Dis_matrix=abs(Dis_matrix);
A=(Dis_matrix<=Cretical_dis);%A:adjacency matrix

figure(1),clf(1),hold on
gplot(A,[Coor_x,Coor_y],'green');
scatter(Coor_x(m+1:end),Coor_y(m+1:end),15,'blue o')
plot(Coor_x(1:m),Coor_y(1:m),'blue*');
axis equal
xlim([0,1])
ylim([0 1])
box on
% 
% gin=ginput;
% ginlen=length(gin(:,1));
% n=n+ginlen;
% Coor_x=[Coor_x;gin(:,1)];
% Coor_y=[Coor_y;gin(:,2)];

clear i;
Coor=Coor_x+i*Coor_y;
Dis_matrix=Coor*ones(1,n)-ones(n,1)*conj(Coor');
Dis_matrix=abs(Dis_matrix);
A=(Dis_matrix<=Cretical_dis);%A:adjacency matrix
Connect=A;
Connect(1:m,:)=0;
Connect=Connect-diag(diag(Connect));
Connect_num=sum(sum(Connect));
figure(7),spy(Connect)

Sensorreading_Dis_matrix=A.*Cretical_dis;
figure(1),clf(1),hold on
gplot(A,[Coor_x,Coor_y],'green');
scatter(Coor_x(m+1,end),Coor_y(m+1,end),15,'blue')
plot(Coor_x(1:m),Coor_y(1:m),'blue*');
axis equal
xlim([0,1])
ylim([0 1])
box on

% dt=.005*10^(-5);
% t0=0;
% ts=10*10^(-5);


dt=.05*10^(-5);
t0=0;
ts=10*10^(-5);


epsilon2=1*10^5;
epsilon1=1*10^5;
f=zeros(n,n);
dfdx1=zeros(n,n);
dfdx2=zeros(n,n);


Lambda0=rands(n,n).*A;
E_Coor_x0=rand(n,1);
E_Coor_y0=rand(n,1);






Lambda=Lambda0;
E_Coor_x=E_Coor_x0;
E_Coor_y=E_Coor_y0;

E_Coor_x(1:m)=Coor_x(1:m);
E_Coor_y(1:m)=Coor_y(1:m);

xstore=E_Coor_x0;
ystore=E_Coor_y0;

error_store=zeros(1,0);
tic
for t=t0:dt:ts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %For a particular problem, this part may differ.
    dfdx1=2*(E_Coor_x*ones(1,n)-ones(n,1)*E_Coor_x');
    dfdx2=2*(E_Coor_y*ones(1,n)-ones(n,1)*E_Coor_y');
    f=(E_Coor_x*ones(1,n)-ones(n,1)*E_Coor_x').^2+(E_Coor_y*ones(1,n)-ones(n,1)*E_Coor_y').^2-Sensorreading_Dis_matrix.^2;
    f_new=f.*Connect;
    f_error=sum(sum(f_new.*(f_new>0)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %MYminus=((f+Lambda)<=0).*(f+Lambda);
    MYminus=((f+Lambda)>=0).*(f+Lambda);
    
    dx1=-epsilon1*dt*diag(Lambda*dfdx1');
    dx2=-epsilon1*dt*diag(Lambda*dfdx2');
    
    dLambda=-epsilon2*dt*(Lambda-MYminus);
    
    E_Coor_x=E_Coor_x+dx1;
    E_Coor_y=E_Coor_y+dx2;
    
    Lambda=Lambda+dLambda;
    Lambda=Lambda.*A;
    
    
E_Coor_x(1:m)=Coor_x(1:m);
E_Coor_y(1:m)=Coor_y(1:m);

xstore=[xstore,E_Coor_x];
ystore=[ystore,E_Coor_y];
error_store=[error_store,f_error];
end
toc
    f=(E_Coor_x*ones(1,n)-ones(n,1)*E_Coor_x').^2+(E_Coor_y*ones(1,n)-ones(n,1)*E_Coor_y').^2-Sensorreading_Dis_matrix.^2;
    f_new=f.*Connect;
    f_error=sum(sum(f_new.*(f_new>0)));
    error_store=[error_store,f_error];


% figure(1),hold on
% plot(E_Coor_x,E_Coor_y,'ro')

figure(2),hold on,
for i=m+1:n
h=plot3(t0:dt:ts+dt,xstore(i,:),ystore(i,:));
%set(h,'linewidth',2);
box on
end
xlabel('time/s');ylabel('x');zlabel('y');
xlim([t0,ts]);

figure(3),clf(3),hold on,
for i=1:n
    plot([E_Coor_x(i),Coor_x(i)],[E_Coor_y(i),Coor_y(i)],'green-')
end
plot(E_Coor_x(m+1:end),E_Coor_y(m+1:end),'ro')
scatter(Coor_x,Coor_y,15,'blue')
plot(Coor_x(1:m),Coor_y(1:m),'blue*');
axis equal
xlim([0,1])
ylim([0 1])
box on

figure(4),clf(4),hold on,
plot(t0:dt:ts+dt,error_store/Connect_num)
xlim([t0,ts])
box on,


figure(5),hold on,
h=plot(t0:dt:ts+dt,xstore(m+1:end,:));
%set(h,'linewidth',2);
box on
ylabel('x')
xlabel('t/s')
xlim([t0,ts])


figure(6),hold on,
h=plot(t0:dt:ts+dt,ystore(m+1:end,:));
%set(h,'linewidth',2);
box on
ylabel('y')
xlabel('t/s')
xlim([t0,ts])


myerror=sqrt(((norm(E_Coor_x(m+1:end)-Coor_x(m+1:end)))^2+(norm(E_Coor_y(m+1:end)-Coor_y(m+1:end))).^2)/(n-m))