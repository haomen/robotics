clear all
close all
display('===========')
%parameters
%sensor net parameters
n=2;%node number
range=0.95;%communication range
eta=0.55;%10 percent ranging error
eta_minus=1-1/(eta+1);
eta_plus=1/(1-eta)-1;
delta=10*2*pi/360;%10 degree bearing error
%algorithm parameters
ts=500;
c0=0.00050;
k0=1;
k1=1;
k2=1;
k3=1;
%node ID to visualize
n_vis=1;
%generate nodes
node_coor=[0.4756 + 0.4554i;   0.5699 + 0.5523i];
%%%%%%%%%%%
dis_matrix=[0    0.1352;    0.1352         0];
angle_matrix=[0    0.7990;    -2.3425         0];
comm_matrix=[0 1;1 0];%this matrix indicates the communication topology
Dis_measure_matrix=[         0    0.1129;    0.1129         0];
Angle_measure_matrix=[0    0.8361;    -2.3055         0];
%Angle_measure_matrix=(Angle_measure_matrix+(pi+Angle_measure_matrix').*comm_matrix)/2;
%%%%%%%%%%%
%visualize the topology
figure(1),clf(1),hold on
gplot(comm_matrix,[real(node_coor),imag(node_coor)],'green');
hold on
scatter(real(node_coor),imag(node_coor),15);
axis equal
xlim([-0.15,1.15])
ylim([-0.15,1.15])
box on
title('WSN topology')
fig_handle=figure(2);clf(2),hold on
plot(real(node_coor(n_vis)),imag(node_coor(n_vis)),'sr');
text(real(node_coor(n_vis)+0.007),imag(node_coor(n_vis)),num2str(n_vis));
axis equal
box on
for i=1:n
    if i==n_vis
        continue;
    elseif comm_matrix(n_vis,i)>0
  figure(2),hold on
plot(real(node_coor(i)),imag(node_coor(i)),'black+');
plot([real(node_coor(n_vis)) real(node_coor(i))],[imag(node_coor(n_vis)) imag(node_coor(i))],'green');
text(real(node_coor(i))+0.007,imag(node_coor(i)),num2str(i));
d0=Dis_measure_matrix(n_vis,i);
theta0=Angle_measure_matrix(n_vis,i);
eta01=d0*eta_plus;
eta02=d0*eta_minus;
draw_constraint([real(node_coor(n_vis)),imag(node_coor(n_vis))],d0,eta01,eta02,theta0,delta,fig_handle)
    end
end
%%%%%%%%%%%
%state variables to run the algorithm
x=rands(n,2);%associated with node
%xxxxxxxxxxx
x=[real(node_coor),imag(node_coor)];
%xxxxxxxxxxx
xstore=zeros(n,2,ts+1);
xstore(:,:,1)=x;
lambda=100*rand(n,n,3);%associated with link--only the elements consistent with the topology 9comm_matrix) count.
for i=1:3
lambda(:,:,i)=(lambda(:,:,i)+lambda(:,:,i)')/2;
lambda(:,:,i)=lambda(:,:,i).*comm_matrix;
end
lambdastore=zeros(n,n,3,ts+1);
lambdastore(:,:,:,1)=lambda;
mu=rand(n,n);%associated with link--only the elements consistent with the topology 9comm_matrix) count.
mu=(mu+mu')/2;
mu=mu.*comm_matrix;
mustore=zeros(n,n,ts+1);
mustore(:,:,1)=mu;
    %compute A,b
    A11=-sin(Angle_measure_matrix+delta);
    A11=A11.*comm_matrix;
    A12=cos(Angle_measure_matrix+delta);
    A12=A12.*comm_matrix;
    A21=sin(Angle_measure_matrix-delta);
    A21=A21.*comm_matrix;
    A22=-cos(Angle_measure_matrix-delta);  
    A22=A22.*comm_matrix;
    A31=-cos(Angle_measure_matrix);
    A31=A31.*comm_matrix;
    A32=-sin(Angle_measure_matrix);
    A32=A32.*comm_matrix;
    ANode1=[A11(1,2),A12(1,2);A21(1,2),A22(1,2);A31(1,2),A32(1,2)];
    ANode2=[A11(2,1),A12(2,1);A21(2,1),A22(2,1);A31(2,1),A32(2,1)];
    b3=Dis_measure_matrix.*(1-eta_plus);
    b3=(b3>=0).*b3*cos(delta);
    b3=b3.*comm_matrix;
    
 a12=[-sin(Angle_measure_matrix(1,2)+delta),cos( Angle_measure_matrix(1,2)+delta);...
      sin(Angle_measure_matrix(1,2)-delta),-cos( Angle_measure_matrix(1,2)-delta);...
     -cos(Angle_measure_matrix(1,2)),-sin(Angle_measure_matrix(1,2))];
 a21=[-sin(Angle_measure_matrix(2,1)+delta),cos( Angle_measure_matrix(2,1)+delta);...
      sin(Angle_measure_matrix(2,1)-delta),-cos( Angle_measure_matrix(2,1)-delta);...
     -cos(Angle_measure_matrix(2,1)),-sin(Angle_measure_matrix(2,1))];%xxxxxxxxxxxxxxxxxxxxxxxx    
B12=[0;0;b3(1,2)];
B21=[0;0;b3(2,1)];
%%

%%


 %     A11=0;
%     A12=0;
%     A21=0;
%     A22=0;
% A31=0;
% A32=0;
% b3=0;
%xxxxxxxxxxxxxxxxxxxxxxxx
t_start= cputime;
%iterations
for t=1:ts
    clear i;
    Dis_vec_matrix=-(x(:,1)+x(:,2)*i)*ones(1,n)+ones(n,1)*conj((x(:,1)+x(:,2)*i)');
    Dis_vec_matrix=Dis_vec_matrix.*comm_matrix;
    current_dis_matrix=abs(Dis_vec_matrix);
    Lmu=diag(mu*ones(n,1))-mu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Aij xj-Aij xi
%the result include three rows
%row 1:Aij^11 (xj(1)-xi(1))+Aij^12(xj(2)-xi(2))
%in matrix form,
%row 1:A11.*real(Dis_vec_matrix)+A12.*imag(Dis_vec_matrix);
Ax_row1=A11.*real(Dis_vec_matrix)+A12.*imag(Dis_vec_matrix);
%row 2:Aij^21 (xj(1)-xi(1))+Aij^22(xj(2)-xi(2))
%in matrix form,
%row 2:A21.*real(Dis_vec_matrix)+A22.*imag(Dis_vec_matrix);
Ax_row2=A21.*real(Dis_vec_matrix)+A22.*imag(Dis_vec_matrix);
%row 3:Aij^31 (xj(1)-xi(1))+Aij^32(xj(2)-xi(2))
%in matrix form,
%row 3:A31.*real(Dis_vec_matrix)+A32.*imag(Dis_vec_matrix);
Ax_row3=A31.*real(Dis_vec_matrix)+A32.*imag(Dis_vec_matrix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %sum_{j\in\N(i)mu_{ij}(x_j-x_i)=-L(mu)x
    Lmu_x=-Lmu*x;%two column vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %sum_{j\in\N(i)A_{ij}^T lambda_{ij}= a two-row vector: 
    %row 1. sum_{j\in\N(i)[A^{11}_ij lambda^1_{ij}+A^{21}_ij lambda^2_{ij}+A^{31}_ij lambda^3_{ij}]
    %In matrix form, 
    %row 1 column vector=trace(A11'*lambda(:,:,1)+A21'*lambda(:,:,2)+A31'*lambda(:,:,3));
    lambda(:,:,1)=lambda(:,:,1).*comm_matrix;
    lambda(:,:,2)=lambda(:,:,2).*comm_matrix;
    lambda(:,:,3)=lambda(:,:,3).*comm_matrix;
    
    A_lambda_row1_column=((A11'*lambda(:,:,1)+A21'*lambda(:,:,2)+A31'*lambda(:,:,3)))*ones(n,1);
    %row 2. sum_{j\in\N(i)[A^{12}_ij lambda^1_{ij}+A^{22}_ij lambda^2_{ij}+A^{32}_ij lambda^3_{ij}]
    %In matrix form,
    %row 2 column vector=trace(A12'*lambda(:,:,1)+A22'*lambda(:,:,2)+A32'*lambda(:,:,3));
    A_lambda_row2_column=((A12'*lambda(:,:,1)+A22'*lambda(:,:,2)+A32'*lambda(:,:,3)).*comm_matrix)*ones(n,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% proj_omega
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%update
x
mu
lambda

xnext=x;
xnext(:,1)=(1-c0)*x(:,1)+c0*proj_omega(-inf,inf,...
    x(:,1)+k1*(k0*A_lambda_row1_column)+k1*Lmu_x(:,1)...
);
xnext(:,2)=(1-c0)*x(:,2)+c0*proj_omega(-inf,inf,...
    x(:,2)+k1*(k0*A_lambda_row2_column)+k1*Lmu_x(:,2)...
);
Ax_row=zeros(n,n,3);
Ax_row(:,:,1)=Ax_row1;
Ax_row(:,:,2)=Ax_row2;
Ax_row(:,:,3)=Ax_row3;
b_matrix=zeros(n,n,3);
b_matrix(:,:,3)=b3;
lambdanext=(1-c0)*lambda+c0*proj_pos(...
    lambda+k2*Ax_row+k2*b_matrix...
);
munext=(1-c0)*mu+c0*proj_pos(...
    mu+k3*current_dis_matrix.^2 ...
    -k3*(Dis_measure_matrix.*(1+eta_plus)).^2);

display('hi')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
x=xnext;
lambda=lambdanext;
mu=munext;
%xxxxxxxxxxxxxxxxxxxxxxxx
%mu=0;
%xxxxxxxxxxxxxxxxxxxxxxxx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%save data
xstore(:,:,t+1)=x;
lambdastore(:,:,:,t+1)=lambda;
mustore(:,:,t+1)=mu;
end
t_time= cputime-t_start;
display(['total time:',num2str(t_time),' seconds'])
figure(3),clf(3),hold on,box on
xdata=0:size(xstore,3)-1;
for i=1:n
xdata(1:end)=xstore(i,1,1:end);
subplot(2,1,1),hold on,
plot(0:size(xstore,3)-1,xdata);
title('x'),box on
xdata(1:end)=xstore(i,2,1:end);
subplot(2,1,2),hold on,
plot(0:size(xstore,3)-1,xdata)
title('y'),box on
end

figure(4),clf(4),hold on,box on
for i=1:n
    for j=1:n
        for k=1:3
           if (comm_matrix(i,j)==0)|(j==i)
               continue
           end
xdata(1:end)=lambdastore(i,j,k,:);
subplot(3,1,k),hold on
plot(0:length(xdata)-1,xdata);title(['lambda',num2str(k)]),box on
        end
    end
end

figure(5),clf(5),hold on,box on
for i=1:n
    for j=1:n
           if (comm_matrix(i,j)==0)|(j==i)
               continue
           end
xdata(1:end)=mustore(i,j,:);
plot(0:length(xdata)-1,xdata);
    end
end
title('mu')
figure(6),clf(6),hold on,box on
xdata=0:size(xstore,3)-1;
xdata1=xdata;
for i=1:n
xdata(1:end)=xstore(i,1,1:end);
xdata1(1:end)=xstore(i,2,1:end);
plot(xdata,xdata1,'b')
plot(xdata(1),xdata1(1),'+')
plot(xdata(end),xdata1(end),'o')
end
axis equal
% xlim([-0.15,1.15])
% ylim([-0.15,1.15])
box on
title('evolution of the estimated position')