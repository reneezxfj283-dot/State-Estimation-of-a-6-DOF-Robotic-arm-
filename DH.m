%% DH坐标系建立

syms q1 q2 q3 q4 q5 q6 real
global T_all
% 定义改进DH参数：[alpha a d theta],每一行对应一个连杆，每一列对应一个DH参数（alpha a d theta）单位：m\rad
% dhparams = [
%  %  a       alpha       d        theta 
%     0         0       0.1675      q1;       % joint1
%     0       pi/2      0.1448      q2;      % joint2
%   0.626       0         0         q3;       % joint3
%   0.5755      0       0.0245      q4;      % joint4
%     0      -pi/2      0.106       q5;       % joint5
%     0       pi/2      0.1026      q6        % joint6
% ];

dhparams = [
 %  a       alpha       d        theta 
    0         0       0.1607      q1;       % joint1
    0       pi/2      0.138       q2;      % joint2
  0.626       0       -0.1237     q3;       % joint3
  0.5755      0       0.099       q4;      % joint4
    0      -pi/2      0.099       q5;       % joint5
    0       pi/2      0.0936      q6        % joint6
];

% 齐次变换矩阵的含变量表达式
mdh_transform = @(alpha,a,d,theta) [...
    cos(theta), -sin(theta), 0, a;
    sin(theta)*cos(alpha), cos(theta)*cos(alpha), -sin(alpha), -d*sin(alpha);
    sin(theta)*sin(alpha), cos(theta)*sin(alpha),  cos(alpha),  d*cos(alpha);
    0, 0, 0, 1];

T_all = sym(zeros(4,4,6));
for i = 1:size(dhparams,1)
    a     = dhparams(i,1);
    alpha = dhparams(i,2);
    d     = dhparams(i,3);
    theta = dhparams(i,4);
    T_all(:,:,i) = simplify(mdh_transform(alpha, a, d, theta));
end

disp('--- 各关节的齐次变换矩阵 (符号表达式) ---');
for i = 1:6
    fprintf('T_%d = \n', i);
    disp(T_all(:,:,i));
end

% 创建机器人模型
robot = robotics.RigidBodyTree('DataFormat','column','MaxNumBodies',6);
% 创建并添加刚体
prevBodyName = 'base';     % 记录前一个刚体名称（连接用），base指初始坐标系名称，根节点（基座）
dhparams_num = double(subs(dhparams, [q1 q2 q3 q4 q5 q6], zeros(1,6)));
for i = 1:size(dhparams,1)     % 循环每个 DH 行
    body = robotics.RigidBody(['link' num2str(i)]);               % 创建一个刚体（连杆）
    joint = robotics.Joint(['joint' num2str(i)], 'revolute');     % 创建此连杆与上一个连杆间的旋转关节

    % 使用改进DH参数直接设置变换
    setFixedTransform(joint, dhparams_num(i,:), 'mdh');      % mdh:表示改进DH变换
    body.Joint = joint;                   % 将关节绑定到刚体（连杆）
    addBody(robot, body, prevBodyName);   % 将刚体添加到机器人模型
    prevBodyName = body.Name;             % 添加完每一个link，将名字存入prevBodyName
end
config = homeConfiguration(robot);
config(2) = pi/2; %+
config(4) = -pi/2; %-
figure(1)
showdetails(robot)
show(robot,config,'Frames','on', 'PreservePlot', false);
view(135, 20);
axis([-0.5,0.5,-0.5,0.5,-0.3,1.7])
% title('六轴机械臂各关节坐标系（改进DH法）');
set(gcf,'Units','centimeters','Position',[5 5 12 9]);% 美化
set(gcf,'Color','w');
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1);
h = findobj(gca,'Type','line');
set(h,'LineWidth',1.5);
exportgraphics(gcf,'sim_result.pdf','ContentType','vector');% 同时存矢量 + 高分辨率 PNG
exportgraphics(gcf,'sim_result.png','Resolution',600);


