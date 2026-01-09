function [x,P,u,Out_x,Out_P]= IEE_algorithm(x,P,yk,Model,u)
TP= Model.trs;               % transition matrix
M = size(TP,1);              % number of model
N = size(Model.Ad{1},1);       % 2   % dimension of state
A = Model.Ad;
B = Model.Bd;
C = Model.H;
Q = Model.Q;
R = Model.R;
w = [200;200;200;55;55;55]; 
mid_x = zeros(N,M);
mid_P = zeros(N,N,M);
for j=1 : M
 % % Mixing
    Pred_u = TP(:,j)'.* u; % TP(:,j)'前一时刻的模型转移概率矩阵第j列转置；u前一时刻的模型权重向量

 % Prediction
    Pred_x = A{j} * x + B{j} * w;   
    Pred_P = A{j} * P * A{j}' + 0.35*eye(12);
    S   = C{j}  * Pred_P * C{j}'+ R{j};
    lihood_u(j)= sum(Pred_u) * mvnpdf(yk,C{j}*Pred_x,S);

 % Approximate
    mid_x(:,j) = C{j}' * R{j}^-1 * yk + Pred_P^-1 * Pred_x;
    mid_P(:,:,j) = C{j}' * R{j}^-1 * C{j} + Pred_P^-1;
end
    
    u=lihood_u/sum(lihood_u);    % update model probability
    inv_P = zeros(N,N);
    for i=1 : M
        summP =  u(i) * mid_P(:,:,i);
        inv_P =  inv_P + summP;
    end
   P = inv_P^-1;
   %P = pinv(inv_P);
   Out_P = P;
   x = Out_P * sum(mid_x.* repmat(u,N,1),2);
   Out_x = x;
  
end