clc
clear
%% creating the W matrix
%V=load('J16.mat').a
%J=load('V16.mat').b
V=load('synthetic_V.mat').random_V
%V=load('v1732.mat').v
J=load('synthetic_J.mat').J
inv_V = inv(V)
V_prim = transpose(V)
inv_V_prim = inv(V_prim)
W = V*J*inv_V
%% build P
muu=zeros(16,1)
D_out = sum(W,2)
D_out = diag(D_out.^-1)
P = D_out * W
P_prim = (P+eye(16))/2
b=[muu ; 1]
a=P_prim' - eye(16)
a_p=[a;ones(1,16)]

x=linsolve(a_p,b)
st_s = diag(x)
L_P = st_s - (st_s*P_prim + (P_prim)'*st_s)/2

%% extracting R_mds and redefining s_S
R_mds = extracted_Rmd(J)
size_rmd = size(R_mds)
%% undirected settings
%W1 = W + W'

%[u1 landa1] = eig(L_P)
s_fun = @(x) exp(-x)
S = filter_definition(R_mds,s_fun)
snr_o=[]
snr_z=[]
snr_t=[]
%% define H
h_fun = @(x) exp(-0.5*x)
H_J = filter_definition(R_mds,h_fun)
H = V*H_J*inv_V


%% start
s_S = S*(S')
new_S_s = s_S
mu = zeros(1,16)
covx = V*V_prim
cov_wn = covx

for snr_in =1:1:1
    snr_o1=[]
    snr_z1=[]
    snr_t1=[]
    for i=1:1:1
        %% signal producing
        w_s = mvnrnd(mu,covx,1)
        signal =  V*S*inv_V *(w_s')

        w_n= mvnrnd(mu,cov_wn,1)
        a_c=10^(snr_in/20)
        d=20*log10(norm(H*signal)/norm(w_n))
        c = 10^((d/20) - (snr_in/20))
        wn = c*w_n'
        sigma = c
        %% noise properties

        G_n = sigma^2* eye(size(J))



        Beta =  5*10^-1
%% determining epsilon for tikhonov which relates to cov of white noise 
        new_cov = sigma^2 * cov_wn     %cause original w_n had cov of V*V' but wn has cov of sigma^2*V*V'
        epsilon = sqrt(trace(new_cov )) 











        %% define y
        y = H*signal + wn
        y_hat = inv_V*y
        %% ینصمص
        beta1 = Beta
        x_l = []
        g_S=[]
        w_J=[]
        g_s_j=[]
        g_J=[]
        F_matrix=[]
        Z_zegond_matrix=[]
        for i =1:length(R_mds)-1
            r_md_i=R_mds(:,i)
            sizee = r_md_i(2) - r_md_i(1) + 1
            Hl = H_J(r_md_i(1):r_md_i(2) , r_md_i(1):r_md_i(2))
            P = (Hl') * (Hl)
            B = new_S_s( r_md_i(1):r_md_i(2) , r_md_i(1):r_md_i(2))
            g_n= G_n(r_md_i(1):r_md_i(2) , r_md_i(1):r_md_i(2))
            C = (Hl')*g_n*(Hl)
            T = P*B*P + C
            E = P*B + B*P
            F = kron(eye(sizee) , T) + kron(transpose(T) , eye(sizee))
            g = reshape(E,[],1)
            vec_z_zegond = inv(F) * g
            z_zegond = reshape(vec_z_zegond,sizee,sizee)
            F_matrix = blkdiag(F_matrix, rank(F) - size(F,1))
            Z_zegond_matrix=blkdiag(Z_zegond_matrix ,rank(z_zegond)-size(z_zegond,1))
            Z = inv(z_zegond) - P   %% W_hermitian * W
            x_l = [x_l ; z_zegond*(Hl')*y_hat(r_md_i(1):r_md_i(2),1)]
            w_J = blkdiag(w_J ,Z)

            g_s_j = inv(beta1*Z +eye(sizee))
            g_J = blkdiag(g_J, g_s_j)
        end
        g_S = V * g_J * inv(V)
        %% recover the signal
        sigma_xy = V*s_S*V_prim*(H')
        sigma_y = H * V*s_S*V_prim*(H') + (sigma^2)*V*V_prim

        x_init = rand(16,1)
        u_init= x_init
        %recovered1 = weiner_optimization(x_init ,u_init, g_S,V,V_prim,inv_V,inv_V_prim,H,Beta,y,signal)
        y_hat = inv_V * y
        recovered2 = V * recovery_mmse(R_mds , y_hat,w_J, H_J)
        recovered3 = recovery_mmse_invertex(H,V,V_prim,inv_V,inv_V_prim,w_J,y)
        %undirected_recovery = weiner_optimization2(x_init,H,g_L,y)
        snr_Zero = 20*log10( norm(inv(H)*y) / norm(inv(H)*y - signal) )
        %snr_out=(norm(signal)/norm(signal - recovered1))
        %l1 = recovered1 - recovered2
        l2 = recovered2 - recovered3
        x_hat = sigma_xy * inv(sigma_y)*y
        snr_out2=20*log10(norm(signal)/norm(signal - x_hat))

        snr_o1=[snr_o1 snr_out2]
        snr_z1=[snr_z1 snr_Zero]

        o= optimvar('o',16);
        prob = optimproblem;
        prob.Objective = transpose(o) * L_P * o
        prob.Constraints.cons1 = norm(o -y)<=epsilon
        x0.o=zeros(1,16)
        sol=solve(prob,x0)
        sol_o=sol.o
        snr_tik = 20*log10(norm(signal)/norm(signal - sol.o))
        snr_t1=[snr_t1 snr_tik]
        %snr_undi = 20*log10(norm(signal)/norm(signal - undirected_recovery))
        %l4 = max(abs(V*V_prim*inv(H')*inv_V_prim*w_J*inv_V-(sigma^2)*V*V_prim*inv(H')*inv_V_prim*inv(new_S_s)*inv_V),[],'all')
    end
    snr_o=[snr_o ; snr_o1]
    snr_t=[snr_t ; snr_t1]
    snr_z=[snr_z ; snr_z1]
end
snr_directed = mean(snr_o,2)
snr_tikhonov=mean(snr_t,2)
snr_zero_forcing=mean(snr_z ,2)
paramplot.position = [100,100,300,220];
plot(ss , snr_tikhonov ,ss,snr_directed ,ss,snr_zero_forcing,'LineWidth',2)
xlabel('Input SNR (db) ');
ylabel('Output SNR (dB)');
axis tight;
title('deconvolution on Synthetic Data')
legend('Directed Tikhonov','Directed Wiener','Zero Forcing','Location','Best');
%% the function that calculates the covariance matrix
function cov_X = cov_x(x,number_of_cities)
l=number_of_cities
cov=zeros(l,l)
for i=1:l
    for j=1:l
        cov(i,j) = sum(x(i,:).*x(j,:),'all') / length(x)

    end
end
cov_X=cov

end

%% extracting blocks
function extracted_Rmd = extracted_Rmd (J)
l=length(J)
x = []
R_md=[]
for i = 1: l
    if i ~= l  &&  J(i,i+1) ==0
        x=[x i]

    end

end
x=[0 x]
for i = 1:length(x)
    if i~=length(x)
        R_md = [R_md [x(i)+1 x(i+1) J(x(i+1),x(i+1))]']
    else
        R_md = [R_md [x(i)+1 l J(l,l)]']
    end
end
extracted_Rmd = R_md
end
%% defining block diagonal H filter based on a function
function H = filter_definition(R_md , func)
%recieves mathematical function and R_mds and returns the block diagonal H
Block_d=[]
len=length(R_md)
for i = 1:len
    r_md_i=R_md(:,i)
    r=r_md_i(2,1)-r_md_i(1,1)
    a=[]
    a(1)=func(r_md_i(3,1))
    if r==0

        Block_d=blkdiag(Block_d,a)
    else
        for j =1:r
            syms x
            g=diff(func(x),j)
            g=matlabFunction( g )
            try
                a(j+1) = g(r_md_i(3,1))./factorial(j)
            catch
                a(j+1)=feval(g)

            end
        end
        h=triu(toeplitz(a))
        Block_d=blkdiag(Block_d,h)
    end

end
H=Block_d
end
%% extracting submatrix
function extracted_submatrix = extracting_submatrix(A,l)
extracted_submatrix = A(l(1):l(2) , l(1):l(2))
end

%% weiner optimization for directed graphs
function recovered_signal = weiner_optimization(x_init , u_init,g_S,V,V_prim,inv_V,inv_V_prim,H,beta,y,signal)
x_j = x_init
u_j=u_init
t_j = 1
for j =1:2000

    v_prim = x_j - beta*V*(V_prim)*(H')*inv_V_prim*inv_V*(H*x_j - y)
    u_j_1=g_S*v_prim
    t_j_1 = (1+sqrt(1+4*(t_j)^2))/2
    x_new = u_j_1 + ((t_j-1)/(t_j_1))*(u_j_1 - u_j)
    x_old = x_j
    x_j = x_new
    u_j = u_j_1
    t_j = t_j_1



end
recovered_signal = x_new
end

%% recovery from 41 relation in paper
function recovery_mmse = recovery_mmse(R_mds, y_hat,Z,H_I)
x_l_hat=[]
for i =1:length(R_mds)
    rmd = R_mds(:,i)
    z = Z(rmd(1,1):rmd(2,1) , rmd(1,1):rmd(2,1))
    h_block = H_I(rmd(1,1):rmd(2,1) , rmd(1,1):rmd(2,1))
    h = (h_block')*h_block
    x_l_recovered = inv(z + h) * (h_block') * y_hat(rmd(1,1):rmd(2,1) , 1)

    x_l_hat = [x_l_hat ; x_l_recovered]

end
recovery_mmse = x_l_hat

end

%% %% recovery from 41 relation in paper in vertex domain

function recovery_mmse_invertex = recovery_mmse_invertex(H,V,V_prim,inv_V,inv_V_prim,Z,y)
recovery_mmse_invertex = inv(H+V*V_prim*inv(H')*inv_V_prim*Z*inv_V)*y
%recovery_mmse_invertex = inv(H+10^-2*(V*V_prim*inv(H')*inv_V_prim*Z*inv_V))*y
end

%% weiner for undirected graph

function weiner_optimization2 = weiner_optimization2 ( x0 , H,  g_L ,y)
% initialization of parameters
z_j=x0
u_j=x0
t_j=[1]


beta = 0.25



for j = 1:2000
    V = z_j - beta* (H')*(H*z_j - y)
    u_j_1 = g_L*V
    t_j_1=(1+sqrt(1+4*t_j^2) )/2
    z_j_1 = u_j_1 + ( (t_j - 1) / (t_j_1) ) *(u_j_1 - u_j)
    z_j =  z_j_1
    u_j = u_j_1
    t_j = t_j_1

end
weiner_optimization2 =  z_j

end