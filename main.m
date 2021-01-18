%clc;
clear all;

% load('yale_mtv.mat');  % gt:165x1, 15 clusters; X: 4096 3304 6750(3)x165
% load('HOG_yale.mat');  % 2700x165
% X{4}=H;
 views = 3;
% load('yaleB_mtv.mat');  % gt:650x1, 10 clusters; X: 2500 3304 6750(3)x165
 load('ORL_mtv.mat');  % gt:400x1, 40 clusters; X: 4096 3304 6750(3)
% load('NH_interval9_mtv.mat');  % gt:550x1, 5 clusters; X: 2000 3304 6750(3)
% load('bbcsport_2view.mat');  % gt:544x1, 5 clusters; X: 3183 3203(2)

iter = 100;
M = 4;
nmix = zeros(M,1);
accr = zeros(M,1);
spa = zeros(M,1);
pur = zeros(M,1);
Fscore = zeros(M,1);
RandIdx = zeros(M,1);
Prec = zeros(M,1);
Rec = zeros(M,1);
ar = zeros(M,1);
err=zeros(iter,1);

lambda_c = 0.05;
lambda_d = 0.1;
%       c      d 
% yale 0.05    1
% orl  0.05  0.1
% nh    0.5    1
% bbc   0.1    1

for m = 1:M
tol=10^-7; 
clusters = size(unique(gt),1);
N = size(gt,1);
mu = 10^-6;
rho=1.1;
maxmu=10^6;
    C = zeros(N,N);
    Y2 = zeros(N,N);
    for v = 1:views
        X{v} = X{v}./repmat(sqrt(sum(X{v}.^2,1)),size(X{v},1),1);
        D{v} = rand(N,N);
        E{v} = zeros(size(X{v}));
        W{v} = zeros(size(X{v}));
        Y1{v} = zeros(size(X{v}));
        Y3{v} = zeros(size(X{v}));
    end
 
     x = zeros(views, 1);
    for i = 1:views
        x(i) = norm(X{i},'fro')^2;
    end
    xx=sum(x);
    for s = 1:iter
        tic;
       
        for i=1:views
            %% D
            Ad = zeros(N,N);
            Cd = zeros(N,N);
            Ad = mu*X{i}'*X{i}+2*lambda_d.*eye(N);
            Cd = mu*X{i}'*(X{i}+Y1{i}/mu-X{i}*C-E{i});
            D{i} = Ad\Cd;
            
            %% E
            E{i} = (X{i}+Y1{i}/mu-X{i}*C-X{i}*D{i}+W{i}-Y3{i}/mu)/2;
            
            %% W
            temp = E{i}+Y3{i}/mu;
            W{i} = solve_l1l2(temp,1/mu);
        end
        
        %% K
        [U,S,V] = svd(C + Y2./mu);
        a = diag(S)-lambda_c/mu;
        a(a<0)=0;
        T = diag(a);
        K = U*T*V';
        
        %% C
        A = zeros(N,N);
        B = zeros(N,N);
        for i=1:views
            A = A + X{i}'*X{i};
            B = B + X{i}'*(X{i}+Y1{i}/mu-X{i}*D{i}-E{i});
        end
        C = (A+eye(N,N))\(B+K-Y2/mu);

        %% check if converge
        for i = 1:views
            p(i) = norm(X{i}-X{i}*(C+D{i}),'fro')^2;
            eq2{i} = X{i}-X{i}*(C+D{i})-E{i};
            eq3{i} = E{i}-W{i};
        end
        P = sum(p);
        eq1=C-K;
        err(s)=abs(xx-P);
        %eq=[norm(eq1,Inf),norm(eq2{1},Inf),norm(eq2{2},Inf),norm(eq2{3},Inf),err(s)];
% bbcsport        
%eq=[norm(eq1,Inf),norm(eq2{1},Inf),norm(eq2{2},Inf),err(s)];
	%stop=max(eq);
        if s==iter || err(s)<tol
            break;
        else
            xx = P;
            for i=1:views
                Y1{i}=Y1{i}+mu*eq2{i};
                Y3{i}=Y3{i}+mu*eq3{i};
            end
            Y2 = Y2 + mu*eq1;
            mu = min(maxmu,mu*rho);
        end
        s = s+1;
    end
    
    %% clustering
    Zt=zeros(N,N);
    for i = 1:views
        Zt=Zt+(abs(D{i})+abs(D{i})');
    end
    Z = (abs(C)+abs(C'))/2+Zt/(2*views);
    Clus = SpectralClustering(Z,clusters);
    time = toc;
    
end  
