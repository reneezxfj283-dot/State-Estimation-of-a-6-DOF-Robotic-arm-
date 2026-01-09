function [x,P,u,Out_x,Out_P]= IMM_algorithm(x,P,y,Model,u)
TP= Model.trs;
M = size(TP,1);
N = size(Model.Ad{1},1);      % 2   % dimension of state
A = Model.Ad;
B = Model.Bd;
C = Model.H;
Q = Model.Q;
R = Model.R;
w=[190;190;190;40;30;40];
upd_x = zeros(N,M);
upd_P = zeros(N,N,M);
for j=1 : M
    
    Pred_u = TP(:,j)'.* u;         
    Norm_u = Pred_u/sum(Pred_u);
    Mix_x  = sum(x.*repmat(Norm_u,N,1),2);     % mixed state vector
    Mix_P  = zeros(N,N);
    for i =1 : M
    summP =  Norm_u(i)*(P(:,:,i)+(x(:,i)-Mix_x)*(x(:,i)-Mix_x)');
    Mix_P =  Mix_P + summP;
    end                                     % mixed covariance
    Pred_x = A{j} * Mix_x + B{j} * w;
    Pred_P = A{j} * Mix_P * A{j}' + 0.35*eye(12);
    S   = C{j}  * Pred_P * C{j}'+ R{j};
    K   = Pred_P * C{j}' * S^-1; 
    lihood_u(j)= sum(Pred_u) * mvnpdf(y,C{j}*Pred_x,S);
    upd_x(:,j) = Pred_x + K *(y - C{j} * Pred_x);
    upd_P(:,:,j)   = Pred_P - K * C{j} * Pred_P;
   
end
    u=lihood_u/sum(lihood_u);
    x=upd_x;
    P=upd_P;
    Out_x = sum(x.* repmat(u,N,1),2);      % Output
    Out_P = zeros(N,N);
    for i=1 : M
    summP =  u(i)*(P(:,:,i)+(x(:,i)-Out_x)*(x(:,i)-Out_x)');
    Out_P =  Out_P + summP;   
    end
    
end