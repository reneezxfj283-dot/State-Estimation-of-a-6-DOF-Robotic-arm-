clc;
clear;
close all

M         = 2;
sim_number = 1000;
steps      = 250;%仿真总步数
rob_Para;

%% Initial values

P0 = diag([10,10,10,10,10,10,100,100,100,100,100,100]);
Init_u=[1/2 1/2];

% Storage Data
Mse_kf1 = zeros(sim_number,steps);    Mse_kf2 = zeros(sim_number,steps);    Mse_kf3 = zeros(sim_number,steps);   
Mse_kf4 = zeros(sim_number,steps);    Mse_kf5 = zeros(sim_number,steps);    Mse_kf6 = zeros(sim_number,steps);  
Mse_kf7 = zeros(sim_number,steps);    Mse_kf8 = zeros(sim_number,steps);    Mse_kf9 = zeros(sim_number,steps);   
Mse_kf10 = zeros(sim_number,steps);   Mse_kf11 = zeros(sim_number,steps);   Mse_kf12 = zeros(sim_number,steps);   

Mse_imm1 = zeros(sim_number,steps);   Mse_imm2 = zeros(sim_number,steps);   Mse_imm3 = zeros(sim_number,steps);  
Mse_imm4 = zeros(sim_number,steps);   Mse_imm5 = zeros(sim_number,steps);   Mse_imm6 = zeros(sim_number,steps);  
Mse_imm7 = zeros(sim_number,steps);   Mse_imm8 = zeros(sim_number,steps);   Mse_imm9 = zeros(sim_number,steps);   
Mse_imm10 = zeros(sim_number,steps);  Mse_imm11 = zeros(sim_number,steps);  Mse_imm12 = zeros(sim_number,steps);    

Mse_iee1 = zeros(sim_number,steps);   Mse_iee2 = zeros(sim_number,steps);   Mse_iee3 = zeros(sim_number,steps);   
Mse_iee4 = zeros(sim_number,steps);   Mse_iee5 = zeros(sim_number,steps);   Mse_iee6 = zeros(sim_number,steps);   
Mse_iee7 = zeros(sim_number,steps);   Mse_iee8 = zeros(sim_number,steps);   Mse_iee9 = zeros(sim_number,steps);   
Mse_iee10 = zeros(sim_number,steps);  Mse_iee11 = zeros(sim_number,steps);  Mse_iee12 = zeros(sim_number,steps);    

for cas = 1  :  1 
Mdl = Rob_Model();

for simN = 1 : sim_number                    % Monto-Carlo Simulation
m0 = [0; 0; 2*pi/3; -pi/4; -pi/2; 0; 0; 0; 0; 0; 0; 0];  % initial state                      
t=0;
x = Gauss(m0,P0);
r=rand(1); 
if r<=Init_u(1)
    m=1;
% elseif Init_u(1)<r&&r<=Init_u(1)+Init_u(2)
%     m=2;
else
    m=2;
end

X = zeros(size(m0,1),steps);
S = zeros(1,steps);
Y = zeros(size(Mdl.H{1},1),steps);
T = zeros(1,steps);

for k=1 : steps
     r=rand(1); 
%  Determin the mode sequences     
    if  m==1
        if r<=Mdl.trs(1,1)
            m=1;
        % elseif Mdl.trs(1,1)<r && r<=Mdl.trs(1,1)+Mdl.trs(1,2)
        %     m=2;
        else
            m=2;
        end
    elseif m==2
        if r<=Mdl.trs(2,1)
            m=1;
        % elseif Mdl.trs(2,1)<r && r<= Mdl.trs(2,1)+Mdl.trs(2,2)
        %     m=2;
        else
            m=2;
        end
    % else
    %     if r<=Mdl.trs(3,1)
    %         m=1;
    %     elseif Mdl.trs(3,1)<r && r<= Mdl.trs(3,1)+Mdl.trs(3,2)
    %         m=2;
    %     else
    %         m=3;
    %     end
    end

    x = Mdl.Ad{m} * x + Gauss(zeros(12,1),Mdl.Q{m})+ Mdl.Bd{m} * Mdl.I{m} ;
    y = Mdl.H{m} * x + Gauss(zeros(6,1),Mdl.R{m});
    
    t = t + Ts; 
    X(:,k) = x;
    Y(:,k) = y;
    S(:,k) = m;
    T(k)   = t;    
end

%%  KF Algoithm
KF_x = zeros(size(m0,1),size(Y,2));
KF_P = zeros(size(m0,1),size(m0,1),size(Y,2));

x = m0; 
P=P0;
 for k = 1 : size(Y,2)  
 yk = Y(:,k);
 [x,P] = KF_algorithm(x,P,yk,Mdl,1);   
 KF_x(:,k)  = x;
 KF_P(:,:,k)= P;  

 % Out_x=m0;
 % Out_P=P0;
 % u = Init_u;
 % x = repmat(m0,1,M);
 % P = ones(size(m0,1),size(m0,1),M); 
 % [x,P,u,Out_x,Out_P] = KF_algorithm(x,P,yk,Mdl,u,Out_x,Out_P);  
 % KF_x(:,k)  = Out_x;
 % KF_P(:,:,k)= Out_P;  
 end
KF_error = (KF_x -X).^2;
Mse_kf1(simN,:) =  KF_error(1,:);    Mse_kf2(simN,:) = KF_error(2,:);     Mse_kf3(simN,:) =  KF_error(3,:);       
Mse_kf4(simN,:) = KF_error(4,:);     Mse_kf5(simN,:) =  KF_error(5,:);    Mse_kf6(simN,:) = KF_error(6,:);
Mse_kf7(simN,:) =  KF_error(7,:);    Mse_kf8(simN,:) = KF_error(8,:);     Mse_kf9(simN,:) =  KF_error(9,:);       
Mse_kf10(simN,:) = KF_error(10,:);   Mse_kf11(simN,:) =  KF_error(11,:);  Mse_kf12(simN,:) = KF_error(12,:);

%% IMM Algorithm
IMM_x = zeros(size(m0,1),size(Y,2));
IMM_P = zeros(size(m0,1),size(m0,1),size(Y,2));
IMM_u = zeros(M,size(Y,2));

u = Init_u;
x = repmat(m0,1,M);
P = ones(size(m0,1),size(m0,1),M);
for i=1:M
     P(:,:,i)= P0 ;
end

 for k = 1 : size(Y,2)   
 [x,P,u,Out_x,Out_P]= IMM_algorithm(x,P,Y(:,k),Mdl,u);   % IMM Algorithm

  IMM_x(:,k)  = Out_x;
  IMM_P(:,:,k)= Out_P;
  IMM_u(:,k) = u';   
 end
IMM_error = (IMM_x -X).^2;
Mse_imm1(simN,:) =  IMM_error(1,:);    Mse_imm2(simN,:) = IMM_error(2,:);      Mse_imm3(simN,:) =  IMM_error(3,:);       
Mse_imm4(simN,:) = IMM_error(4,:);     Mse_imm5(simN,:) =  IMM_error(5,:);     Mse_imm6(simN,:) = IMM_error(6,:);
Mse_imm7(simN,:) =  IMM_error(7,:);    Mse_imm8(simN,:) = IMM_error(8,:);      Mse_imm9(simN,:) =  IMM_error(9,:);       
Mse_imm10(simN,:) = IMM_error(10,:);   Mse_imm11(simN,:) =  IMM_error(11,:);   Mse_imm12(simN,:) = IMM_error(12,:);

%% IEE Algorithm
IEE_x = zeros(size(m0,1),size(Y,2));
IEE_P = zeros(size(m0,1),size(m0,1),size(Y,2));
IEE_u = zeros(M,size(Y,2));

u = Init_u;
x=m0;
P=P0;

 for k = 1 : size(Y,2) 
 yk = Y(:,k);
 [x,P,u,Out_x,Out_P]= IEE_algorithm(x,P,yk,Mdl,u); 

  IEE_x(:,k)  = Out_x;
  IEE_P(:,:,k)= Out_P;
  IEE_u(:,k) = u';   
 end
IEE_error = (IEE_x -X).^2;
Mse_iee1(simN,:) =  IEE_error(1,:);   Mse_iee2(simN,:) = IEE_error(2,:);      Mse_iee3(simN,:) =  IEE_error(3,:);       
Mse_iee4(simN,:) = IEE_error(4,:);    Mse_iee5(simN,:) =  IEE_error(5,:);     Mse_iee6(simN,:) = IEE_error(6,:);
Mse_iee7(simN,:) =  IEE_error(7,:);   Mse_iee8(simN,:) = IEE_error(8,:);      Mse_iee9(simN,:) =  IEE_error(9,:);       
Mse_iee10(simN,:) = IEE_error(10,:);  Mse_iee11(simN,:) =  IEE_error(11,:);   Mse_iee12(simN,:) = IEE_error(12,:);

end
end

%% RMSE_Velocity
figure(1)
subplot(2,2,1)
plot(1:steps,sqrt(mean(Mse_iee8)),'m',1:steps,sqrt(mean(Mse_imm8)),'b',1:steps,sqrt(mean(Mse_kf8)),'g');
legend('Proposed method','IMM','KF','Orientation','horizontal');
legend('boxoff')
% grid on
xlabel('Time(k)');ylabel('joint2 RMSEs');
set(gcf,'Units','centimeters','Position',[5 5 12 9]);  % Enhance color
set(gcf,'Color','w');
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
h = findobj(gca,'Type','line');
set(h,'LineWidth',1.5);
exportgraphics(gcf,'sim_result.pdf','ContentType','vector');
exportgraphics(gcf,'sim_result.png','Resolution',600);

subplot(2,2,2)
plot(1:steps,sqrt(mean(Mse_iee9)),'m',1:steps,sqrt(mean(Mse_imm9)),'b',1:steps,sqrt(mean(Mse_kf9)),'g');
legend('Proposed method','IMM','KF','Orientation','horizontal');
legend('boxoff')
% grid on
xlabel('Time(k)');ylabel('joint3 RMSEs');
set(gcf,'Units','centimeters','Position',[5 5 12 9]);% Enhance color
set(gcf,'Color','w');
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
h = findobj(gca,'Type','line');
set(h,'LineWidth',1.5);
exportgraphics(gcf,'sim_result.pdf','ContentType','vector');
exportgraphics(gcf,'sim_result.png','Resolution',600);

subplot(2,2,3)
plot(1:steps,sqrt(mean(Mse_iee11)),'m',1:steps,sqrt(mean(Mse_imm11)),'b',1:steps,sqrt(mean(Mse_kf11)),'g');
legend('Proposed method','IMM','KF','Orientation','horizontal');
legend('boxoff')
% grid on
xlabel('Time(k)');ylabel('joint5 RMSEs');
set(gcf,'Units','centimeters','Position',[5 5 12 9]);% Enhance color
set(gcf,'Color','w');
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
h = findobj(gca,'Type','line');
set(h,'LineWidth',1.5);
exportgraphics(gcf,'sim_result.pdf','ContentType','vector');
exportgraphics(gcf,'sim_result.png','Resolution',600);

subplot(2,2,4)
plot(1:steps,sqrt(mean(Mse_iee12)),'m',1:steps,sqrt(mean(Mse_imm12)),'b',1:steps,sqrt(mean(Mse_kf12)),'g');
legend('Proposed method','IMM','KF','Orientation','horizontal');
legend('boxoff')
% grid on
xlabel('Time(k)');ylabel('joint6 RMSEs');
set(gcf,'Units','centimeters','Position',[5 5 12 9]);% Enhance color
set(gcf,'Color','w');
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
h = findobj(gca,'Type','line');
set(h,'LineWidth',1.5);
exportgraphics(gcf,'sim_result.pdf','ContentType','vector');
exportgraphics(gcf,'sim_result.png','Resolution',600);

%% RMSE _Position
% figure(1)
% subplot(2,2,1)
% plot(1:steps,sqrt(mean(Mse_iee2)),'m',1:steps,sqrt(mean(Mse_imm2)),'b',1:steps,sqrt(mean(Mse_kf2)),'g');
% legend('Proposed method','IMM','KF','Orientation','horizontal');
% legend('boxoff')
% % grid on
% xlabel('Time(k)');ylabel('joint2 RMSEs');
% set(gcf,'Units','centimeters','Position',[5 5 12 9]);% Enhance color
% set(gcf,'Color','w');
% set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
% h = findobj(gca,'Type','line');
% set(h,'LineWidth',1.5);
% exportgraphics(gcf,'sim_result.pdf','ContentType','vector');
% exportgraphics(gcf,'sim_result.png','Resolution',600);
% 
% subplot(2,2,2)
% plot(1:steps,sqrt(mean(Mse_iee3)),'m',1:steps,sqrt(mean(Mse_imm3)),'b',1:steps,sqrt(mean(Mse_kf3)),'g');
% legend('Proposed method','IMM','KF','Orientation','horizontal');
% legend('boxoff')
% % grid on
% xlabel('Time(k)');ylabel('joint3 RMSEs');
% set(gcf,'Units','centimeters','Position',[5 5 12 9]);% Enhance color
% set(gcf,'Color','w');
% set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
% h = findobj(gca,'Type','line');
% set(h,'LineWidth',1.5);
% exportgraphics(gcf,'sim_result.pdf','ContentType','vector');
% exportgraphics(gcf,'sim_result.png','Resolution',600);
% 
% subplot(2,2,3)
% plot(1:steps,sqrt(mean(Mse_iee4)),'m',1:steps,sqrt(mean(Mse_imm4)),'b',1:steps,sqrt(mean(Mse_kf4)),'g');
% legend('Proposed method','IMM','KF','Orientation','horizontal');
% legend('boxoff')
% % grid on
% xlabel('Time(k)');ylabel('joint4 RMSEs');
% set(gcf,'Units','centimeters','Position',[5 5 12 9]);% Enhance color
% set(gcf,'Color','w');
% set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
% h = findobj(gca,'Type','line');
% set(h,'LineWidth',1.5);
% exportgraphics(gcf,'sim_result.pdf','ContentType','vector');
% exportgraphics(gcf,'sim_result.png','Resolution',600);
% 
% subplot(2,2,4)
% plot(1:steps,sqrt(mean(Mse_iee5)),'m',1:steps,sqrt(mean(Mse_imm5)),'b',1:steps,sqrt(mean(Mse_kf5)),'g');
% legend('Proposed method','IMM','KF','Orientation','horizontal');
% legend('boxoff')
% % grid on
% xlabel('Time(k)');ylabel('joint5 RMSEs');
% set(gcf,'Units','centimeters','Position',[5 5 12 9]);% Enhance color
% set(gcf,'Color','w');
% set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
% h = findobj(gca,'Type','line');
% set(h,'LineWidth',1.5);
% exportgraphics(gcf,'sim_result.pdf','ContentType','vector');
% exportgraphics(gcf,'sim_result.png','Resolution',600);

%% ARMSE_Velocity
% RMSE_KF7 = sqrt(mean(Mse_kf7));     RMSE_KF8 = sqrt(mean(Mse_kf8));      RMSE_KF9 = sqrt(mean(Mse_kf9));
% RMSE_KF10 = sqrt(mean(Mse_kf10));   RMSE_KF11 = sqrt(mean(Mse_kf11));    RMSE_KF12 = sqrt(mean(Mse_kf12));
% RMSE_IMM7 = sqrt(mean(Mse_imm7));   RMSE_IMM8 = sqrt(mean(Mse_imm8));    RMSE_IMM9 = sqrt(mean(Mse_imm9));
% RMSE_IMM10 = sqrt(mean(Mse_imm10)); RMSE_IMM11 = sqrt(mean(Mse_imm11));  RMSE_IMM12 = sqrt(mean(Mse_imm12));
% RMSE_IEE7 = sqrt(mean(Mse_iee7));   RMSE_IEE8 = sqrt(mean(Mse_iee8));    RMSE_IEE9 = sqrt(mean(Mse_iee9));     
% RMSE_IEE10 = sqrt(mean(Mse_iee10)); RMSE_IEE11 = sqrt(mean(Mse_iee11));  RMSE_IEE12 = sqrt(mean(Mse_iee12));
% 
% % 逐步计算 ARMSE
% for t = 1:steps
%     ARMSE_KF8(t) = mean(RMSE_KF8(1:t));  ARMSE_KF10(t) = mean(RMSE_KF10(1:t));  ARMSE_KF12(t) = mean(RMSE_KF12(1:t));
%     ARMSE_KF7(t) = mean(RMSE_KF7(1:t));  ARMSE_KF9(t) = mean(RMSE_KF9(1:t));  ARMSE_KF11(t) = mean(RMSE_KF11(1:t));
%     ARMSE_IMM8(t) = mean(RMSE_IMM8(1:t));  ARMSE_IMM10(t) = mean(RMSE_IMM10(1:t));  ARMSE_IMM12(t) = mean(RMSE_IMM12(1:t)); 
%     ARMSE_IMM7(t) = mean(RMSE_IMM7(1:t));  ARMSE_IMM9(t) = mean(RMSE_IMM9(1:t));  ARMSE_IMM11(t) = mean(RMSE_IMM11(1:t));
%     ARMSE_IEE8(t) = mean(RMSE_IEE8(1:t));  ARMSE_IEE10(t) = mean(RMSE_IEE10(1:t));  ARMSE_IEE12(t) = mean(RMSE_IEE12(1:t)); 
%     ARMSE_IEE7(t) = mean(RMSE_IEE7(1:t));  ARMSE_IEE9(t) = mean(RMSE_IEE9(1:t));  ARMSE_IEE11(t) = mean(RMSE_IEE11(1:t)); 
% end
% steps = [75, 150, 225, 300];
% 
% ARMSE_7 = [ARMSE_KF7(75), ARMSE_KF7(150), ARMSE_KF7(225), ARMSE_KF7(300);ARMSE_IMM7(75), ARMSE_IMM7(150), ARMSE_IMM7(225), ARMSE_IMM7(300);ARMSE_IEE7(75), ARMSE_IEE7(150), ARMSE_IEE7(225), ARMSE_IEE7(300)]';
% ARMSE_8 = [ARMSE_KF8(75), ARMSE_KF8(150), ARMSE_KF8(225), ARMSE_KF8(300);ARMSE_IMM8(75), ARMSE_IMM8(150), ARMSE_IMM8(225), ARMSE_IMM8(300);ARMSE_IEE8(75), ARMSE_IEE8(150), ARMSE_IEE8(225), ARMSE_IEE8(300)]';
% ARMSE_9 = [ARMSE_KF9(75), ARMSE_KF9(150), ARMSE_KF9(225), ARMSE_KF9(300);ARMSE_IMM9(75), ARMSE_IMM9(150), ARMSE_IMM9(225), ARMSE_IMM9(300);ARMSE_IEE9(75), ARMSE_IEE9(150), ARMSE_IEE9(225), ARMSE_IEE9(300)]';
% ARMSE_10 = [ARMSE_KF10(75), ARMSE_KF10(150), ARMSE_KF10(225), ARMSE_KF10(300);ARMSE_IMM10(75), ARMSE_IMM10(150), ARMSE_IMM10(225), ARMSE_IMM10(300);ARMSE_IEE10(75), ARMSE_IEE10(150), ARMSE_IEE10(225), ARMSE_IEE10(300)]';
% ARMSE_11 = [ARMSE_KF11(75), ARMSE_KF11(150), ARMSE_KF11(225), ARMSE_KF11(300);ARMSE_IMM11(75), ARMSE_IMM11(150), ARMSE_IMM11(225), ARMSE_IMM11(300);ARMSE_IEE11(75), ARMSE_IEE11(150), ARMSE_IEE11(225), ARMSE_IEE11(300)]';
% ARMSE_12 = [ARMSE_KF12(75), ARMSE_KF12(150), ARMSE_KF12(225), ARMSE_KF12(300);ARMSE_IMM12(75), ARMSE_IMM12(150), ARMSE_IMM12(225), ARMSE_IMM12(300);ARMSE_IEE12(75), ARMSE_IEE12(150), ARMSE_IEE12(225), ARMSE_IEE12(300)]';
% figure(1);
% subplot(2,3,1)
% bar(steps, ARMSE_7, 'grouped');
% xlabel('Steps');  % x轴标签
% ylabel('ARMSE of joint1');  % y轴标签
% legend('KF', 'IMM', 'Proposed method');ylim([0 20]);
% grid on;  % 设置网格
% subplot(2,3,2)
% bar(steps, ARMSE_8, 'grouped');
% xlabel('Steps');  % x轴标签
% ylabel('ARMSE of joint2');  % y轴标签
% legend('KF', 'IMM', 'Proposed method');ylim([0 60]);
% grid on;  % 设置网格
% subplot(2,3,3)
% bar(steps, ARMSE_9, 'grouped');
% xlabel('Steps');  % x轴标签
% ylabel('ARMSE of joint3');  % y轴标签
% legend('KF', 'IMM', 'Proposed method');ylim([0 60]);
% grid on;  % 设置网格
% subplot(2,3,4)
% bar(steps, ARMSE_10, 'grouped');
% xlabel('Steps');  % x轴标签
% ylabel('ARMSE of joint4');  % y轴标签
% legend('KF', 'IMM', 'Proposed method');ylim([0 30]);
% grid on;  % 设置网格
% subplot(2,3,5)
% bar(steps, ARMSE_11, 'grouped');
% xlabel('Steps');  % x轴标签
% ylabel('ARMSE of joint5');  % y轴标签
% legend('KF', 'IMM', 'Proposed method');ylim([0 400]);
% grid on;  % 设置网格
% subplot(2,3,6)
% bar(steps, ARMSE_12, 'grouped');
% xlabel('Steps');  % x轴标签
% ylabel('ARMSE of joint6');  % y轴标签
% legend('KF', 'IMM', 'Proposed method');ylim([0 200]);
% grid on;  % 设置网格

%% ARMSE_Position
% RMSE_KF2 = sqrt(mean(Mse_kf2));  RMSE_KF4 = sqrt(mean(Mse_kf4));  RMSE_KF6 = sqrt(mean(Mse_kf6));
% RMSE_IMM2 = sqrt(mean(Mse_imm2));  RMSE_IMM4 = sqrt(mean(Mse_imm4));  RMSE_IMM6 = sqrt(mean(Mse_imm6));
% RMSE_IEE2 = sqrt(mean(Mse_iee2));  RMSE_IEE4 = sqrt(mean(Mse_iee4));  RMSE_IEE6 = sqrt(mean(Mse_iee6));
% RMSE_KF1 = sqrt(mean(Mse_kf1));  RMSE_KF3 = sqrt(mean(Mse_kf3));  RMSE_KF5 = sqrt(mean(Mse_kf5));
% RMSE_IMM1 = sqrt(mean(Mse_imm1));  RMSE_IMM3 = sqrt(mean(Mse_imm3));  RMSE_IMM5 = sqrt(mean(Mse_imm5));
% RMSE_IEE1 = sqrt(mean(Mse_iee1));  RMSE_IEE3 = sqrt(mean(Mse_iee3));  RMSE_IEE5 = sqrt(mean(Mse_iee5));
% for t = 1:steps
%     ARMSE_KF2(t) = mean(RMSE_KF2(1:t));  ARMSE_KF4(t) = mean(RMSE_KF4(1:t));  ARMSE_KF6(t) = mean(RMSE_KF6(1:t));
%     ARMSE_IMM2(t) = mean(RMSE_IMM2(1:t));  ARMSE_IMM4(t) = mean(RMSE_IMM4(1:t));  ARMSE_IMM6(t) = mean(RMSE_IMM6(1:t));  
%     ARMSE_IEE2(t) = mean(RMSE_IEE2(1:t));  ARMSE_IEE4(t) = mean(RMSE_IEE4(1:t));  ARMSE_IEE6(t) = mean(RMSE_IEE6(1:t)); 
%     ARMSE_KF1(t) = mean(RMSE_KF1(1:t));  ARMSE_KF3(t) = mean(RMSE_KF3(1:t));  ARMSE_KF5(t) = mean(RMSE_KF5(1:t));
%     ARMSE_IMM1(t) = mean(RMSE_IMM1(1:t));  ARMSE_IMM3(t) = mean(RMSE_IMM3(1:t));  ARMSE_IMM5(t) = mean(RMSE_IMM5(1:t)); 
%     ARMSE_IEE1(t) = mean(RMSE_IEE1(1:t));  ARMSE_IEE3(t) = mean(RMSE_IEE3(1:t));  ARMSE_IEE5(t) = mean(RMSE_IEE5(1:t));      
% end
% steps = [75, 150, 225, 300];
% ARMSE_1 = [ARMSE_KF1(75), ARMSE_KF1(150), ARMSE_KF1(225), ARMSE_KF1(300);
%     ARMSE_IMM1(75), ARMSE_IMM1(150), ARMSE_IMM1(225), ARMSE_IMM1(300);
%     ARMSE_IEE1(75), ARMSE_IEE1(150), ARMSE_IEE1(225), ARMSE_IEE1(300)]';
% ARMSE_2 = [ARMSE_KF2(75), ARMSE_KF2(150), ARMSE_KF2(225), ARMSE_KF2(300);
%     ARMSE_IMM2(75), ARMSE_IMM2(150), ARMSE_IMM2(225), ARMSE_IMM2(300);
%     ARMSE_IEE2(75), ARMSE_IEE2(150), ARMSE_IEE2(225), ARMSE_IEE2(300)]';
% ARMSE_3 = [ARMSE_KF3(75), ARMSE_KF3(150), ARMSE_KF3(225), ARMSE_KF3(300);
%     ARMSE_IMM3(75), ARMSE_IMM3(150), ARMSE_IMM3(225), ARMSE_IMM3(300);
%     ARMSE_IEE3(75), ARMSE_IEE3(150), ARMSE_IEE3(225), ARMSE_IEE3(300)]';
% ARMSE_4 = [ARMSE_KF4(75), ARMSE_KF4(150), ARMSE_KF4(225), ARMSE_KF4(300);
%     ARMSE_IMM4(75), ARMSE_IMM4(250), ARMSE_IMM4(225), ARMSE_IMM4(300);
%     ARMSE_IEE4(75), ARMSE_IEE4(150), ARMSE_IEE4(225), ARMSE_IEE4(300)]';
% ARMSE_5 = [ARMSE_KF5(75), ARMSE_KF5(150), ARMSE_KF5(225), ARMSE_KF5(300);
%     ARMSE_IMM5(75), ARMSE_IMM5(150), ARMSE_IMM5(225), ARMSE_IMM5(300);
%     ARMSE_IEE5(75), ARMSE_IEE5(150), ARMSE_IEE5(225), ARMSE_IEE5(300)]';
% ARMSE_6 = [ARMSE_KF6(75), ARMSE_KF6(150), ARMSE_KF6(225), ARMSE_KF6(300);
%     ARMSE_IMM6(75), ARMSE_IMM6(150), ARMSE_IMM6(225), ARMSE_IMM6(300);
%     ARMSE_IEE6(75), ARMSE_IEE6(150), ARMSE_IEE6(225), ARMSE_IEE6(300)]';
% figure(1);
% subplot(1,3,1)
% bar(steps, ARMSE_4, 'grouped');
% xlabel('Steps');  % x轴标签
% ylabel('ARMSE of joint 4');  % y轴标签
% legend('KF', 'IMM', 'Proposed method');
% grid on;  % 设置网格
% subplot(1,3,2)
% bar(steps, ARMSE_5, 'grouped');
% xlabel('Steps');  % x轴标签
% ylabel('ARMSE of joint 5');  % y轴标签
% legend('KF', 'IMM', 'Proposed method');
% grid on;  % 设置网格
% subplot(1,3,3)
% bar(steps, ARMSE_6, 'grouped');
% xlabel('Steps');  % x轴标签
% ylabel('ARMSE of joint 6');  % y轴标签
% legend('KF', 'IMM', 'Proposed method');
% grid on;  % 设置网格

%% true-filter plot
% figure(1);
% subplot(1,2,1);
% plot(1:steps, X(6,:), '--', 'LineWidth', 1.3); hold on;
% plot(1:steps, KF_x(6,:),'g', 'LineWidth', 1.3);hold on;
% plot(1:steps, IMM_x(6,:),'b', 'LineWidth', 1.3);hold on;
% plot(1:steps, IEE_x(6,:),'m', 'LineWidth', 1.3);
% legend('true','KF','IMM','Proposed method','Orientation','horizontal');
% xlabel('time(k)');
% ylabel('Angle');
% set(gcf,'Units','centimeters','Position',[5 5 12 9]);% 美化
% set(gcf,'Color','w');
% set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
% h = findobj(gca,'Type','line');
% set(h,'LineWidth',1.5);
% exportgraphics(gcf,'sim_result.pdf','ContentType','vector');% 同时存矢量 + 高分辨率 PNG
% exportgraphics(gcf,'sim_result.png','Resolution',600);
% 
% subplot(1,2,2);
% plot(1:steps, X(12,:), '--', 'LineWidth', 1.3); hold on;
% plot(1:steps, KF_x(12,:),'g', 'LineWidth', 1.3);hold on;
% plot(1:steps, IMM_x(12,:),'b', 'LineWidth', 1.3);hold on;
% plot(1:steps, IEE_x(12,:),'m', 'LineWidth', 1.3);
% legend('true','KF','IMM','Proposed method','Orientation','horizontal');
% xlabel('time(k)');
% ylabel('Angular Velocity');
% set(gcf,'Units','centimeters','Position',[5 5 12 9]);% 美化
% set(gcf,'Color','w');
% set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
% h = findobj(gca,'Type','line');
% set(h,'LineWidth',1.5);
% exportgraphics(gcf,'sim_result.pdf','ContentType','vector');% 同时存矢量 + 高分辨率 PNG
% exportgraphics(gcf,'sim_result.png','Resolution',600);
 
