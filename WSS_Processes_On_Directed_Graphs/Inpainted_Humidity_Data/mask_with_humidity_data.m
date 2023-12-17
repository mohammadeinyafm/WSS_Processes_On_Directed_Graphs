clc
clear
%% creating Mask filter

N = 10 %number of nodes (cities 1 to 10 between 22 cities)

values = load("h_data.mat").Y % values = humidity data - mean(humidity data)
%% calculate the coviariance of data
cov = cov_x(values,size(values,1)) 
cov_data = cov(1:N,1:N) %cause we need the cov of first 10 cities

%% building graph
J=load('value_2_J.mat').b
V=load('value_2_v.mat').a
inv_V = inv(V)
V_prim = V'
inv_V_prim = inv(V_prim)
W = V*J*inv_V %adjacency of directed matrix
%% build Lp based on the paper that introduced a symetric Laplacian for directed graphs
W_P=W
W_P(1,7)=mean(W,"all")
W_P(3,7)=mean(W,"all")+0.234
D_out = sum(W_P,2)
D_out = diag(D_out.^-1)
P = D_out *W_P
P_prim = (P+eye(N))/2
mu_i=zeros(N,1)
b=[mu_i ; 1]
a=P_prim' - eye(N)
a_p=[a;ones(1,N)]
steady_state=linsolve(a_p,b)
st_s = diag(steady_state)
L_P = st_s - (st_s*P_prim + (P_prim)'*st_s)/2
[U2 landa2]=eig(L_P)
%% rmds (blocks of J with number of columns,rows and also the eigenvalue of that)
R_mds = extracted_Rmd(J)
%% build undirected graph
W=V*J*inv_V
W1=(W+W')/2


[U1 landa] = eig(W1)
%% build S*(S)' based on covariance of data, new_S_s is block shape of s_S (because we have to cut the elements which are not in the blocks)
new_S_s = []
size_rmd = size(R_mds)
s_S = inv_V * cov_data * inv_V_prim
for i=1:size_rmd(1,2)
    extracted_s = extracting_submatrix(s_S , R_mds(:,i))
    new_S_s = blkdiag(new_S_s , extracted_s)
end
%% build S*(S)' for undirected graph (we only hold diagonal of it)
s_S_undi = diag(diag(inv(U1)*cov_data*U1)) % this is for undirected graph
s_S_LP = diag(diag(inv(U2)*cov_data*U2))   % this is for symmetric laplacian Lp


%% covariance of white noise
cov_wn = V*V_prim
mu = zeros(1,N)
%% defining filter H
H1 = eye(N)
H = eye(N)
%% 
snr_o=[]
snr_u=[]
snr_p=[]
snr_t=[]
snr_in=14  
%% start
for i=10:10:90  % the percent of measurement that we want to save
percent =i
    snr_o1=[]
    snr_p1=[]
    snr_u1=[]
    snr_t1=[]

%% start algorithm
    for j=1:10000
%% building a mask filter I_TV based on the percent that we want 
    ind = randperm(N)
    
    Mask = zeros(N,1);
    Mask(ind(1:round(percent/100*N)))=1;
    I_TV = zeros(nnz(Mask),N)
    f=1
    for i=1:N
        if Mask(i,1)==1
            I_TV(f,i) =1
            f = f+1

        end
    end
        random_column = round(unifrnd( 0.5 , 744.5 )) % choosing a random column 
        signal = values(:,random_column)        % choosing a random signal based on the random column
        signal = signal(1:N,1)  %we only need the first 10 data (data on cities of 11 to 32 is not needed)

        w_n= mvnrnd(mu,cov_wn,1) %generate a noise
%% determine a multiplier for w_n in order to have desirable snr_i 
        a_c=10^(snr_in/20)
        d=20*log10(norm(signal)/norm(w_n))
        c = 10^((d/20) - (snr_in/20))

        wn = c*w_n'  
       
        sigma = c % original w_n had the covariance of V*V'-but wn has the covariance of c^2 *V*V'-this c is playing the role of sigma
%% determining epsilon for tikhonov which relates to cov of white noise 
        new_cov = sigma^2 * cov_wn     %cause original w_n had cov of V*V' but wn has cov of sigma^2*V*V'
        epsilon = sqrt(trace(I_TV * new_cov * transpose(I_TV)))
%% G_n is g_n*g_n' in the paper which is equal to sigma^2 * eye
        G_n = (sigma^2)* eye(size(J))  %spectrum of white noise that we generated

        G_n_undi = diag(diag(inv(U1)*new_cov*(U1)))  %spectrum of white noise that we generated but for undirected graph 
        G_n_Lp   =  diag(diag(inv(U2)*new_cov*(U2)))  %spectrum of white noise that we generated but for undirected graph

        %% recovery relations for directed graph
        V_TV = I_TV*V
        U1_TV = I_TV*U1
        U2_TV = I_TV*U2
        P = (H1') * pinv(V_TV) *(V_TV)* H1
        B = new_S_s
        C = (H1') * pinv(V_TV) *(V_TV) *(G_n)* pinv(V_TV) *(V_TV)* (H1)
        E = P*B+B*P
        T= P*B*P+C
        F = kron(eye(N) , T) + kron(transpose(T) , eye(N))
        g = reshape(E,[],1)
        vec_z_zegond = pinv(F) * g
        z_zegond = reshape(vec_z_zegond,N,N)
        Z_prim = pinv(z_zegond)
        Z = Z_prim-P

       
        %% recovery relations for the undirected graph whcih is a special case of directed graph

        P_undi = (H1') * pinv(U1_TV) *(U1_TV)* H1
        B_undi = s_S_undi
        C_undi = (H1') * pinv(U1_TV) *(U1_TV) *(G_n_undi)* pinv(U1_TV) *(U1_TV)* (H1)
        E_undi = P_undi*B_undi+B_undi*P_undi
        T_undi= P_undi*B_undi*P_undi+C_undi
        F_undi = kron(eye(N) , T_undi) + kron(transpose(T_undi) , eye(N))
        g_undi = reshape(E_undi,[],1)
        vec_z_zegond_undi = pinv(F_undi) * g_undi
        z_zegond_undi = reshape(vec_z_zegond_undi,N,N)
        Z_prim_undi = pinv(z_zegond_undi)
        Z_undi = Z_prim-P_undi


 %% recovery relations for Lp symmetric laplacian of directed graph 
        P_LP = (H1') * pinv(U2_TV) *(U2_TV)* H1
        B_LP = s_S_LP
        C_LP = (H1') * pinv(U2_TV) *(U2_TV) *(G_n_Lp)* pinv(U2_TV) *(U2_TV)* (H1)
        E_LP = P_LP*B_LP+B_LP*P_LP
        T_LP= P_LP*B_LP*P_LP+C_LP
        F_LP = kron(eye(N) , T_LP) + kron(transpose(T_LP) , eye(N))
        g_LP = reshape(E_LP,[],1)
        vec_z_zegond_LP = pinv(F_LP) * g_LP
        z_zegond_LP = reshape(vec_z_zegond_LP,N,N)
        Z_prim_LP = pinv(z_zegond_LP)
        Z_LP = Z_prim_LP-P_LP


        %% defining y
        y = H*signal + wn

        y_TV = I_TV * y

%% calculate recovered signal
        recovered = V*z_zegond*(H1')*pinv(V_TV)*y_TV  % recovery for directed

        recovered_undi = U1*z_zegond_undi*(H1')*pinv(U1_TV)*y_TV %recovery for undirected

        recovered_LP = U2*z_zegond_LP*(H1')*pinv(U2_TV)*y_TV %recovery for symmetric Lp

        %% finding solution for tikhonov problem
        o= optimvar('o',N);
        prob = optimproblem;
        prob.Objective = transpose(o) * L_P * o
        prob.Constraints.cons1 = norm(I_TV*o -y_TV)<=epsilon
        x0.o=zeros(1,N)
        sol=solve(prob,x0)
        sol_o=sol.o
        snr_tik = norm(signal)/norm(signal - sol.o)
        snr_t1(end+1)=snr_tik
%% check that every thing is right according to different mathematical relations in the paper. if it is needed you can commend this part
        recovered2 = V*pinv( (H1') * pinv(V_TV) * (V_TV) * (H1) + (Z) ) * (H1')*pinv(V_TV)*y_TV
        sigma_x = V*new_S_s*V_prim
        sigma_xy = sigma_x*((I_TV*H)')
        sigma_y = (I_TV*H)*sigma_x*((I_TV*H)')+(sigma^2)*V_TV*(V_TV)'
        recovered3 = sigma_xy*inv(sigma_y)*y_TV
        snr_out2=(norm(signal)/norm(signal - recovered3))
        l=recovered - recovered2
        l2=recovered3- recovered
%% save the snr_s (not in DB yet)
        snr_out=norm(signal)/norm(signal - recovered)
        snr_undi =norm(signal)/norm(signal - recovered_undi)
        snr_LP =norm(signal)/norm(signal - recovered_LP)
        if 20*log10(snr_out) < 5
            disp(20*log10(snr_out))
        end
        %% 
        snr_o1(end+1) =snr_out
        snr_u1(end+1)=snr_undi
        snr_p1(end+1)=snr_LP

    end
%% every row of matrices in below relates to aquaired snrs(not in DB) based on the percent.for example first row is for 10 percent-second row for 20 pecent and...
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

figure(1)
ss=[10:10:90]
 paramplot.position = [100,100,300,220];
plot(ss , snr_tikhonov ,ss,snr_undirected_recovery,ss,snr_directed_from_Lp,ss,snr_directed_recovery,'LineWidth',2)
xlabel('saved percentage ');
ylabel('Output SNR (dB)');
axis tight;
title('inpainting on real data')
legend('Directed Tikhonov','Undirected Wiener','Wiener On Directed Lp','Directed Wiener','Location','Best');

figure(2)
paramplot.position = [100,100,600,220];
subplot(121)
imagesc(abs(cov_data))
colorbar
title('Covariance matrix');

subplot(122)
a = -10;
disp = 20*log10(abs(s_S));
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





%% extarcing Rmds function
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
%% extracting submatrices

function extracted_submatrix = extracting_submatrix(A,l)
extracted_submatrix = A(l(1):l(2) , l(1):l(2))
end



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