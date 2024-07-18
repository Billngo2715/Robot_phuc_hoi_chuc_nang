
function q_do=IK_do(goc_deg,EndPoints,plot_q,ID_Run)
% INPUT:
% canhTay: chieu dai canh tay
% cangTay: chieu dai cang tay
% OUTPUT:
% q_do: gia tri goc sau khi giai bai toan IK (rad)
% INPUT:
% goc: bo q do duoc(rad)
% EndPoints: toa do diem cuoi do duoc (mm)
%%cac bien chieu dai canh tay, cang tay, cong thuc toa do diem cuoi
global canhTay cangTay rE
%%chuyen don vi degree --> radian
goc=goc_deg*pi/180;

accurat=0.1;
fprintf('The software is running IK \n');

mi=min(goc);
ma=max(goc);
if ID_Run==4 || ID_Run==5 || ID_Run==6
    q_max=repmat([   ma   0    0],1); %gioi han max min cua dong tac
    q_min=repmat( [  mi   0   0],1);
elseif ID_Run==1 || ID_Run==2 || ID_Run==3
    q_max=[pi,pi,pi,pi,pi,pi,pi];
    q_min=[-pi,-pi,-pi,-pi,-pi,-pi,-pi];
end


E= EndPoints ;
n=size(EndPoints,1);% so diem can tinh toan
q_opt=zeros(n,7);% khoi tao ma tran nghiem
q0i=[goc(1,:) 0 0];% khoi tao q0i cho diem dau tien
distance_Error=zeros(n,1);%khoi tao ma tran sai so
generation=zeros(n,1);%ma tran so luong the he
t=zeros(n,1);%ma tran thoi gia tinh toan
tm=0;
%%tinh toan toa do diem cuoi
%%Tinh ma tran DH
RB=1 %do
if RB ==1
    canhTay=340;
    cangTay=270;
    syms a1 a2 a3 a4 a5
    T01=DH(a1,125.5,0,90);
    T12=DH(-a2+pi/2,0,0,-90);
    T23=DH(a3,canhTay,0,-90);
    T34=DH(a4,0,0,90);
    T45=DH(a5+pi/2,cangTay,0,90);
elseif RB == 2
    canhTay=286.9;
    cangTay=336;
    syms a1 a2 a3 a4 a5
    T01=DH(a1,143,0,90);
    T12=DH(-a2+pi/2,7,0,-90);
    T23=DH(a3,canhTay,0,-90);
    T34=DH(0,0,0,90);
    T45=DH(a5+pi/2,cangTay,0,90);
    end
Te=T01*T12*T23*T34*T45;


%%cong thuc tong quat toa do diem cuoi
rE=[Te(1,4), Te(2,4), Te(3,4)];

for e=1:n
    e
    tic %bat dau diem thoi gian
    feval_best = inf;
    k=1;
    while ( (feval_best >= accurat))&& k<10% chay thuat toan dien khi dat duoc do chinh xac
        switch ID_Run
            case 1
                [gbest, fitbest,iter] = PSO_DMP1(q_min,q_max,q0i,EndPoints(e,:),e,accurat,ID_Run);
            case 2
                [gbest, fitbest, iter] = DE_DMP2(q_min,q_max,q0i,EndPoints(e,:),e,accurat,ID_Run);
            case 3
                [gbest, fitbest, iter] = ISADE_arm(1,q_min,q_max,q0i,EndPoints(e,:),e,accurat,x,ID_Run); %1,3,7,12,(13)
            case 4 %Pro PSO
                [gbest, fitbest,iter] = PSO_DMP5dot(q_min,q_max,q0i,EndPoints(e,:),e,accurat,x,ID_Run);
            case 5 %Pro DE
                [gbest, fitbest, iter] = DE_DMP2(q_min,q_max,q0i,EndPoints(e,:),e,accurat,ID_Run);
            case 6 %Pro ISADE
                [gbest, fitbest, iter] = ISADE(1,q_min,q_max,q0i,EndPoints(e,:),e,accurat,x,ID_Run);
            case 7
                [gbest, fitbest, iter] = OBL_ISADE_DMP2(q_min,q_max,q0i,EndPoints(e,:),e,accurat,x);             
        end    
        if fitbest < feval_best
            Xbest = gbest;
            feval_best = fitbest
        end 
        k=k+1;
    end
 
    generation(e)=(iter-1);
    q_opt(e,:)= Xbest;
    e;
    feval_best;
    Xbest;
    if e==1
    q0i=[goc(1,:) 0 0]
    else
    q0i=Xbest;
    end
        
    te(e)=toc;
    t=te;
    tm=tm+te(e);
end

tm
q=q_opt;
q_do=q(:,1:5);
q_do(1,:)=goc(1,:);
if plot_q==1
t1=linspace(0,5,length(goc));
t2=linspace(0,5,length(EndPoints));
fprintf('average exectution time:')
figure(1)
plot(t2,q_opt(:,1),'r-','LineWidth',2);
hold on
plot(t1,goc(:,1),'b-','LineWidth',2);
legend('q1 predict','q2 measure')
hold off

figure(2)
plot(t2,q_opt(:,2)*180/pi,'r-','LineWidth',2);
hold on
plot(t1,goc(:,2)*180/pi,'b-','LineWidth',2);
legend('q2 predict','q2 measure')
hold off

figure(3)
plot(t2,q_opt(:,3),'r-','LineWidth',2);
hold on
plot(t1,goc(:,3),'b-','LineWidth',2);
legend('q3 predict','q3 measure')
hold off

figure(4)
plot(t2,q_opt(:,4)*180/pi,'r-','LineWidth',2);
hold on
plot(t1,goc(:,4)*180/pi,'b-','LineWidth',2);
legend('q4 predict','q4 measure')
hold off

figure(5)
plot(t2,q_opt(:,5)*180/pi,'r-','LineWidth',2);
hold on
plot(t1,goc(:,5)*180/pi,'b-','LineWidth',2);
legend('q5 predict','q5 measure')
hold off



figure(6)

plot3(EndPoints(:,1),EndPoints(:,2),EndPoints(:,3),'b.','LineWidth',2);
xlabel('X(mm)');
ylabel('Y(mm)');
zlabel('Z(mm)');
grid on
hold on 
for j=1:n   
q_do=q_opt(:,1:5);
q1=q(j,1);q2=q(j,2);q3=q(j,3);q4=q(j,4);q5=q(j,5);%;q6=q(j,6);q7=q(j,7);

% toa do diem cuoi theo q_do
x5=double(subs(Te(1,4),[a1 a2 a3 a4 a5],[q1 q2 q3 q4 q5]));
y5=double(subs(Te(2,4),[a1 a2 a3 a4 a5],[q1 q2 q3 q4 q5]));
z5=double(subs(Te(3,4),[a1 a2 a3 a4 a5],[q1 q2 q3 q4 q5]));

R5(j,:)=[x5 y5 z5];

xE=E(j,1);
yE=E(j,2);
zE=E(j,3);
distance_Error(j)=sqrt((x5-xE)^2+(y5-yE)^2+(z5-zE)^2);
plot3(x5,y5,z5,'r.','LineWidth',2);
title('Stage 2 of drinking water task (lifting stage)')

legend('Refernces Endpoints','Simualtion EndPoints');
hold on
 end

% plot3(x5,y5,z5,'r.','LineWidth',2);
 legend('Refernces Endpoints','Simualtion EndPoints');
hold on
%  end
%  hold off


figure(7)
plot(distance_Error,'b','LineWidth',2);
% title('Brushing teeth task')
title('Stage 2 of drinking water task (lifting stage)') %(reaching stage)
% title(' Playing ball task')
% title('Catching the ball task')
legend('Distance Error')
xlabel('Points');
ylabel(' Distance Error (mm)');
hold off

figure(8)
plot(generation,'r','LineWidth',2);
legend('generation')

figure(9)
plot(t,'k','LineWidth',2);
legend('Time')
hold off

fprintf('The software stoped running PSO for DMP \n');
tTB=sum(t)/n
errorTB=sum(distance_Error)/n
generationTb=sum(generation)/n
STD=std(distance_Error)
end




