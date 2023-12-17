clc
clear
load meteo_molene_u.mat;
N=10
N2=4
A = zeros(22,22)
W=zeros(N,N)
 x = info{4};
y = info{3};
z = info{5};
value2 = value - mean(value,2)
for i = 1 : 22
for j = i+1 : 22
   A(i,j) = norm(value2(i,:)-value2(j,:))
end
end

A = A+transpose(A)
A = A(1:N,1:N)
for i =1:N
   
         ls= A(i,:)
        [~, indices] = sort(ls(:), 'ascend');
        ind = indices(1:N2);
        for j=1:length(ind)
            W(i,ind(j,1)) = 1
        end


   
end
dist=A.^2
W = W-diag(diag(W))


Adjacency_matrix = zeros(N,N)
for i=1:N
    for j=1:N
        if z(i) <=z(j) & i~=j & W(i,j)==1
            Adjacency_matrix(i,j)=dist(i,j)


        end
    end
end

syms sigma
%sigma =0.00000056474462026881546086408031573533
fun = @(x) exp(-sigma*x)
l = fun(Adjacency_matrix)
l(l==1)=0
d_out = sum(l,2)
d_in = sum(l,1)
 d_av = (d_out+d_in')/2
 d_ave = sum(d_av,"all")/N
eq = d_ave==4
sigma2 = solve(eq)
fun2 =@(x) exp(-sigma2*x)
l2 = fun2(Adjacency_matrix)
l2=double(l2)
l2(l2==1)=0
d_out = sum(l2,2)
d_in = sum(l2,1)
 d_av = (d_out+d_in')/2
 d_ave = sum(d_av,"all")/N
[v j] = jordan(l2)
