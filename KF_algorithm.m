function [x,P]= KF_algorithm(x,P,yk,Model,mode)
A = Model.Ad{mode};
B = Model.Bd{mode};
C = Model.H{mode};
Q = Model.Q{mode};
R = Model.R{mode};
w = [200;200;200;55;55;55]; 
  %Prediction
    Pred_x = A * x + B * w;                    % Kalman filter
    Pred_P = A* P * A' + 0.35*eye(12);

  %Update
    S = C * Pred_P * C'+ R;
    K = Pred_P * C' * S^-1; 
    upd_x = Pred_x + K *(yk - C * Pred_x);
    upd_P = Pred_P - K * C * Pred_P;

    x=upd_x;
    P=upd_P;
end


% function [x,P,u,Out_x,Out_P]= KF_algorithm(x,P,y,Model,u,Out_x,Out_P)
% TP= Model.trs;
% M = size(TP,1);
% N = size(Model.Ad{1},1);      % 2   % dimension of state
% A = Model.Ad;
% B = Model.Bd;
% C = Model.H;
% Q = Model.Q;
% R = Model.R;
% % Q=[4^2 4^2 4^2];%%difference:note that Q{} is replaced by Q()
% upd_x = zeros(N,M);
% upd_P = zeros(N,N,M);
% for j=1 : M
% 
%     Pred_u = TP(:,j)'.* u;         
%     Norm_u = Pred_u/sum(Pred_u);
%     Mix_x  = sum(x.*repmat(Norm_u,N,1),2);     % mixed state vector
%     Mix_P  = zeros(N,N);
%     for i =1 : M
%     summP =  Norm_u(i)*(P(:,:,i)+(x(:,i)-Mix_x)*(x(:,i)-Mix_x)');
%     Mix_P =  Mix_P + summP;
%     end                                     % mixed covariance
%     Pred_x = A{j} * Mix_x + B{j} * [190;190;190;40;30;40];
%     % Pred_x = A{j} * Mix_x + B{j} * [200;200;200;55;55;55];
%     Pred_P = A{j} * Mix_P * A{j}' + 0.35*eye(12);
%     S   = C{j}  * Pred_P * C{j}'+ R{j};
%     K   = Pred_P * C{j}' * S^-1; 
%     lihood_u(j)= sum(Pred_u) * mvnpdf(y,C{j}*Pred_x,S);
%     upd_x(:,j) = Pred_x + K *(y - C{j} * Pred_x);
%     upd_P(:,:,j)   = Pred_P - K * C{j} * Pred_P;
% 
% end
%     u=lihood_u/sum(lihood_u);
%     if u(1) > u(2)
%     j = 1;
%     else
%     j = 2;
%     end    
%     Pred_x = A{j} * Out_x + B{j} *[200;200;200;55;55;55];                    % Kalman filter
%     Pred_P = A{j}* Out_P * A{j}' + 0.35*eye(12);
%     S = C{j} * Pred_P * C{j}'+ R{j};
%     K = Pred_P * C{j}' * S^-1; 
%     updkf_x = Pred_x + K *(y - C{j} * Pred_x);
%     updkf_P = Pred_P - K * C{j} * Pred_P;
% 
%     Out_x=updkf_x;
%     Out_P=updkf_P;
% 
% 
%     x=upd_x;
%     P=upd_P;
% 
% end