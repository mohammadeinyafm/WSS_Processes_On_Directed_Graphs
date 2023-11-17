clc
clear
A = zeros(32,32)
W=zeros(16,16)
tempts=load('tempt_values').temt_values
for i = 1 : 31
for j = i+1 : 32
   A(i,j) = norm(tempts(i,:)-tempts(j,:))
end
end

A = A+transpose(A)
A=A(17:32,17:32)
for i =1:16
   
        v = A(i,:)
        [~, indices] = sort(v(:), 'ascend');
        ind = indices(1:10);
        for j=1:length(ind)
            W(i,ind(j,1)) = 1
        end


   
end
W = W-diag(diag(W))

data = load('names_and_lat_longs')
tbl= data.z
lats = tbl.lat
longs=tbl.long
dist=A.^2
Adjacency_matrix = zeros(16,16)
for i=1:16
    for j=1:16
        if lats(i+16) <=lats(j+16) & i~=j & W(i,j)==1
            Adjacency_matrix(i,j)=dist(i,j)


        end
    end
end

sigma = 0.00008553860215052712650053462473806
W_D = exp(-sigma.* Adjacency_matrix)
W_D(W_D==1) = 0
%{
new_W = W_D(17:32,17:32)
new_W(2,16) = sum(W_D(17:32,17:32),'all')/(16*16)
mu=zeros(16,1)
D_out = sum(new_W,2)
D_out = diag(D_out.^-1)
P = D_out * new_W
P_prim = (P+eye(16))/2
b=[mu ; 1]
a=P_prim' - eye(16)
a_p=[a;ones(1,16)]

x=linsolve(a_p,b)   
%}