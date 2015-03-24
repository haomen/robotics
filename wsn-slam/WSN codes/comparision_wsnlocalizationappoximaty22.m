% localization in wsn via a distributed dual neural network
clear all
close all
m=9;%number of anchors
%geometry graph
%cretical distance;
Cretical_dis=0.2;
nn=8;
x1=1/nn:1/nn:(1-1/nn);
x2=1/nn:1/nn:(1-1/nn);
[X1,X2]=meshgrid(x1,x2);
% Coor_x=rand(n,1);
% Coor_y=rand(n,1);
[len1,len2]=size(X1);
n=len1*len2;
n=m+n;
Coor_x=zeros(n,1);
Coor_y=zeros(n,1);
Coor_x(m+1:end)=X1(:)+0.02*rands(n-m, 1);
Coor_y(m+1:end)=X2(:)+0.02*rands(n-m, 1);
%number of total sensors

%Coor_x(1:m)=round(Coor_x(1:m));
%Coor_y(1:m)=round(Coor_y(1:m));
Coor_x(1:9)=[1 1 0 0 0.5 0     0.5   1    0.5]';
Coor_y(1:9)=[0 1 1 0 0.5 0.5   1     0.5   0]';

clear i;
Coor=Coor_x+i*Coor_y;
Dis_matrix=Coor*ones(1,n)-ones(n,1)*conj(Coor');
Dis_matrix=abs(Dis_matrix);
A=(Dis_matrix<=Cretical_dis);%A:adjacency matrix

figure(1),clf(1),hold on
gplot(A,[Coor_x,Coor_y],'green');
hold on
scatter(Coor_x(m+1:end),Coor_y(m+1:end),15,'blue o')
hold on
plot(Coor_x(1:m),Coor_y(1:m),'blue*');
axis equal
xlim([0,1])
ylim([0 1])
box on



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

Large_num=1000000;
Sensorreading_Dis_matrix=(A.*Cretical_dis+(A==0)*Large_num);
 
% 
dt=1*10^(-5);
t0=0;
ts=.0020;
% 
% dt=1*10^(-6);
% t0=0;
% ts=.0002;
% 
% dt=.05*10^(-5);
% t0=0;
% ts=2*10^(-5);


epsilon2=10^(4);
epsilon1=10^(4);
f=zeros(n,n);
dfdx1=zeros(n,n);
dfdx2=zeros(n,n);


Lambda0=rand(n,n).*A;
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
mycount=0;
tic
cvx_begin sdp
variable my_x(n-m)
variable my_y(n-m)
E_Coor_x=[Coor_x(1:9);my_x];
E_Coor_y=[Coor_y(1:9);my_y];
Ecombine=[E_Coor_x,E_Coor_y];

m21=length(E_Coor_x);
for kk=1:m21
    for kkk=1:m21
        if kk==kkk
            continue;
        end
    if Sensorreading_Dis_matrix(kk,kkk)>1000
        continue;
    else        
    [eye(2)*Sensorreading_Dis_matrix(kk,kkk),(Ecombine(kk,:)-Ecombine(kkk,:))';(Ecombine(kk,:)-Ecombine(kkk,:)),Sensorreading_Dis_matrix(kk,kkk)]>=0;
    end
    end
end
cvx_end
toc

E_Coor_x=[Coor_x(1:9);my_x];
E_Coor_y=[Coor_y(1:9);my_y];


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



myerror=sqrt(((norm(E_Coor_x(m+1:end)-Coor_x(m+1:end)))^2+(norm(E_Coor_y(m+1:end)-Coor_y(m+1:end))).^2)/(n-m))


figure(3),clf(3),hold on,
plot(E_Coor_x(m+1:end),E_Coor_y(m+1:end),'ro')
scatter(Coor_x,Coor_y,15,'blue')
plot(Coor_x(1:m),Coor_y(1:m),'blue*');
axis equal
xlim([0,1])
ylim([0 1])
box on
legend('position estimation','true position of blind nodes','beacon nodes')
for i=1:n
    plot([E_Coor_x(i),Coor_x(i)],[E_Coor_y(i),Coor_y(i)],'green-')
end


figure(4),clf(4),hold on,
plot(t0:dt:ts+dt,error_store/Connect_num)
xlim([t0,ts])



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