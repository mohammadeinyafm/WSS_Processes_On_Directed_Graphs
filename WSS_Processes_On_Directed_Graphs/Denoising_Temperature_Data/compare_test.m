clc
clear

%% creating the W matrix - we summon the V-J which relates to adjacency directed matrix of cities 17 to 32
lintimeday = load('line_time_day.mat').lintimeday
V=load('jordan_k10.mat').a
J=load('jordan_J_k10.mat').b
inv_V = inv(V)
V_prim = transpose(V)
inv_V_prim = inv(V_prim)
W = V*J*inv_V
%% building Lp which is the symmetric laplacian for directed graphs (based on the paper that i send before)
W_P=W
W_P(2,16)=mean(W,"all")
D_out = sum(W_P,2)
D_out = diag(D_out.^-1)
P = D_out *W_P
P_prim = (P+eye(16))/2
mu_i=zeros(16,1)
b=[mu_i ; 1]
a=P_prim' - eye(16)
a_p=[a;ones(1,16)]
steady_state=linsolve(a_p,b)
st_s = diag(steady_state)
L_P = st_s - (st_s*P_prim + (P_prim)'*st_s)/2
[u1 landa1]=eig(L_P)
%% Building Undirected graph (W+transpose(W)./2)

W1 = (W + W')/2
[u2 landa2] = eig(W1)  %diagonal decomposition of W1 (undirected adjacency matrix)



%%
snr_o=[]  % snr of directed recovery is going to be saved on this
snr_u=[]  % snr of undirected recovery ....
snr_p=[]  %snr of Lp recovery ...
snr_t=[]  % snr of tikhonov.....
%% load temperature data
tempts=load('tempt_values').temt_values  %this is temperature data
%% we make the mean of data to be equal to zero
meani = mean(tempts,2)
for i =1:size(tempts,1)
    temts(i,:) = (tempts(i,:) - meani (i,1))
end
%% determine the covariance of data
cov = cov_x(temts,size(tempts,1))
cov_data = cov(17:32,17:32)  %we only need the cov of cities 17 to 32

%% extracting R_mds (address of blocks in jordan with assosiated eigen value)
R_mds = extracted_Rmd(J)
size_rmd = size(R_mds)
%% cov of white noise and mean of white noise
mu = zeros(1,16)
cov_wn = V*V_prim
%% s_S = S * hermition(S) based on cov.first we build build s_S then we delet elements which are not in the blocks of jordan.we save this in new_S_s
s_S = inv_V * cov_data * inv_V_prim
new_S_s = []
for i=1:size_rmd(1,2)
    extracted_s = extracting_submatrix(s_S , R_mds(:,i))
    new_S_s = blkdiag(new_S_s , extracted_s)
end
%% determine H filter which is eye in this case cause we are gonna compare between directed and undirected and etc.
H_J = eye(16)
H = eye(16)
%% start
for snr =1:1:20 %for different snr_in which is 20*log10(norm(signal)/norm(wn))
    snr_o1=[]
    snr_p1=[]
    snr_u1=[]
    snr_t1=[]
    for i =1:10000
        %% choose signal
random_column = round(unifrnd( 1 , 744 ))
       signal = temts(:,random_column)
        signal = signal(17:32,1)
        %% here we gonna first choose a noise.then we determine a multiplier for noise in order to reach the desired input snr. c is this multiplier
        w_n= mvnrnd(mu,cov_wn,1)
        a_c=10^(snr/20)
        d=20*log10(norm(signal)/norm(w_n))
        c = 10^((d/20) - (snr/20))
        %% since now the added noise is c*w_n, the cov of noise is c^2*V*V'. so c playes the role of sigma
        wn = c*w_n'
        sigma = c

        new_cov = sigma^2 * cov_wn
        %% determining epsilon for tikhonov which relates to trace of cov matrix of noise
        epsilon = sqrt(trace( new_cov))

        %% noise properties . G_n is g_n*g_n' in the paper which is the power spec of white noise (in directed case it is equal to sigma^2*eye

        G_n = (sigma^2)* eye(size(J)) %power spectrum of noise (on directed)
        G_n_undi = diag(diag(inv(u2)*new_cov*(u2))) %power spectrum of noise for undirected setting
        G_n_Lp   =  diag(diag(inv(u1)*new_cov*(u1)))  %spectrum of white noise that we generated but for undirected graph
        %% this beta is needed for recovery algorithms in directed case
        Beta =  0.5

        %% determine the s*s' for undirected case (we remove the elements which are not on the diagonal of it)

        s_L = inv(u2) * cov_data * u2
       diag_sL = diag(s_L)
       g_landa = diag(diag_sL./(diag_sL + 0.25*(diag(G_n_undi))))
        g_L = u2*g_landa*inv(u2)

        %% we do the same thing for Lp
        s_L_2 = inv(u1) * cov_data * u1
       diag_sL_2 = diag(s_L_2)
        g_landa_2 = diag(diag_sL_2./(diag_sL_2 + 0.25*(diag(G_n_Lp))))
        g_L_2 = u1*g_landa_2*inv(u1)





        %% define filter H
        H_J = eye(16)
        H = eye(16)


        %% here for every block we calculate (w_j*w_j') then we save that in W_J matrix

        beta1 = Beta
        x_l = []
        g_S=[]
        w_J=[]
        g_s_j=[]
        g_J=[]
        F_matrix=[]
        Z_zegond_matrix=[]
        for i =1:size(R_mds,2) 
        
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
            Z = inv(z_zegond) - P  
            w_J = blkdiag(w_J ,Z)

            g_s_j = inv(beta1*Z +eye(sizee)) 
            g_J = blkdiag(g_J, g_s_j) %g_J is needed for wiener algorithm ( Dr.Iraji paper page 7-algorithm 2)
        end

        g_S = V * g_J * inv(V) %g_S is needed for wiener algorithm ( Dr.Iraji paper page 7-algorithm 2)


      

        %% define y
        y = H*signal + wn  %wn = c*w_n
       
%% recover signal in directed case wiener( Dr.Iraji paper page 7-algorithm 2)
        x_init = rand(16,1)
        u_init= x_init
      recovered = weiner_optimization(x_init ,u_init, g_S,V,V_prim,inv_V,inv_V_prim,H,Beta,y,signal)  
%% recover signal.undirected case wiener(Vandergheynst algorithm)
      recovered_undi = weiner_optimization2(x_init,H,g_L,y)
%% recover signal from  Vandergheynst algorithm but with Lp laplacian
      recovered_LP = weiner_optimization2(x_init,H,g_L_2,y)
%% recover signal with tikhonov

        o= optimvar('o',16);
        prob = optimproblem;
        prob.Objective = transpose(o) * L_P * o
        prob.Constraints.cons1 = norm(o -y)<=epsilon
        x0.o=zeros(1,16)
        sol=solve(prob,x0)
        sol_o=sol.o

%% determine the aquired snr (not in DB yet)
      snr_out=norm(signal)/norm(signal - recovered) %directed

       snr_undi =norm(signal)/norm(signal - recovered_undi) %undirected
       snr_LP =norm(signal)/norm(signal - recovered_LP) %LP
     snr_tik = norm(signal)/norm(signal - sol.o) %tikhonov

       snr_t1(end+1)=snr_tik
       snr_o1(end+1) =snr_out
       %snr_u1(end+1)=snr_undi
       snr_p1(end+1)=snr_LP

 
%% در نتیجه در صورت نیاز و تست کردن کامنتشو بردارین
        y_hat = inv_V * y
%% recover based on relation 41 in the Dr.Iraji paper (both in vertex and fourier)
        recovered2 = V * recovery_mmse(R_mds , y_hat,w_J, H_J)  
        recovered3 = recovery_mmse_invertex(H,V,V_prim,inv_V,inv_V_prim,w_J,y,sigma)
        l1 = recovered - recovered2
        l2 = recovered2 - recovered3
 %% recover from relation 38  in the Dr.Iraji paper
        sigma_xy = V*new_S_s*V_prim*(H')
        sigma_y = H * V*new_S_s*V_prim*(H') + (sigma^2)*V*V_prim
        x_hat = sigma_xy * inv(sigma_y)*y

       

    end
%% every row of matrices in below relates to aquaired snrs(not in DB) based on the input snr.for example first row is for snr_in = 1 -second row for snr_in=2  and...
snr_o=[snr_o ; snr_o1]
snr_u=[snr_u ; snr_u1]
snr_p=[snr_p;  snr_p1]
snr_t=[snr_t ; snr_t1]

end
%% average the snrs and turn them into DB
snr_directed_recovery=(20*log10(mean(snr_o,2)))'
snr_undirected_recovery=(20*log10(mean(snr_u,2)))'
snr_directed_from_Lp=(20*log10(mean(snr_p,2)))'
snr_tikhonov=(20*log10(mean(snr_t,2)))'
%% plotting the outputs

ss=[1:1:20]
paramplot.position = [100,100,300,220];
plot(ss , snr_tikhonov ,ss,snr_undirected_recovery,ss,snr_directed_from_Lp,ss,snr_directed_recovery,'LineWidth',2)
xlabel('Input SNR (db) ');
ylabel('Output SNR (dB)');
axis tight;
%title('deconvolution on real data')
%legend('Tikhonov on Lp','undirected wiener','directed Lp wiener','directed wiener','Location','Best');
figure(2)
paramplot.position = [100,100,600,220];
subplot(121)
imagesc(abs(cov_data))
colorbar
title('Covariance matrix');

subplot(122)
a = -10;
disp = 10*log10(abs(s_S));
disp(disp<a) = a;
imagesc(disp)
colorbar
imagesc(disp)
colorbar
title('Covariance matrix in Fourier (dB)');


figure(3)
paramplot.position = [100,100,600,220];

imagesc(abs(W))
colorbar
title('Graph Directed weighted adjacency matrix ');

figure(4)
paramplot.position = [100,100,600,220];

imagesc(abs(W1))
colorbar
title('Graph Undirected weighted adjacency matrix ');

figure(5)
plot(lintimeday,tempts(32,:)-273)
title('Temperature over time for NOIRMOUTIER-EN');
axis tight
xlabel('Days')
ylabel('Temperature in degree C')



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

%% extracting submatrix
function extracted_submatrix = extracting_submatrix(A,l)
extracted_submatrix = A(l(1):l(2) , l(1):l(2))
end

%% weiner optimization for directed graphs
function recovered_signal = weiner_optimization(x_init , u_init,g_S,V,V_prim,inv_V,inv_V_prim,H,beta,y,signal)
x_j = x_init
u_j=u_init
t_j = 1
for j =1:250

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
for i =1:size(R_mds,2)
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

function recovery_mmse_invertex = recovery_mmse_invertex(H,V,V_prim,inv_V,inv_V_prim,Z,y,sigma)
%recovery_mmse_invertex = inv(H+V*V_prim*inv(H')*inv_V_prim*Z*inv_V)*y
recovery_mmse_invertex = inv(H+V*V_prim*inv(H')*inv_V_prim*Z*inv_V)*y
end

%% weiner for undirected graph (vandergeinst)

function weiner_optimization2 = weiner_optimization2 ( x0 , H,  g_L ,y)
% initialization of parameters
z_j=x0
u_j=x0
t_j=[1]


beta = 0.25



for j = 1:250
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
