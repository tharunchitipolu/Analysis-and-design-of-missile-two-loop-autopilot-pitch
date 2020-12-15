function RouthApprox=Routh_Approximation(G,r)

% Computes the r-th order Routh Approximation of a given n-th order
% transfer function G, with 1<=r<=n.

  
if ~isa(G,'tf') || ~isscalar(G)
    error('Input needs to be a SISO transfer function');
end
[num,den]=tfdata(G,'v');
D_fact=num(1)/den(1);
num=num-D_fact*den;

num1=num(end:-1:1)/den(1);
den1=den(end:-1:1)/den(1);
n=length(den1)-1;

if ~isreal(r) || (fix(r)~=r) || (r<1) || (r>n)
    error('Invalid value of reduced model order')
end


%Routh Approximation
if mod(n,2)
    A=[den1(1:2:end);den1(2:2:end)];
    B=[num1(1:2:end);num1(2:2:end)];
else
    A=[den1(1:2:end);den1(2:2:end) 0];
    B=[num1(1:2:end);num1(2:2:end) 0];
end
gam(r)=0;del=gam;
gam(1)=A(1,1)/A(2,1);
if gam(1)<=0
    disp('System Unstable. Routh Approximation does not exist');
    RouthApprox=0;
    return
end
for i=3:r+1

    for j=1:(size(A,2)-1)
        A(i,j)=A(i-2,j+1)-gam(i-2)*A(i-1,j+1);
    end
    gam(i-1)=A(i-1,1)/A(i,1);
    if gam(i-1)<=0
        disp('System Unstable. Routh Approximation does not exist');
        RouthApprox=0;
        return
    end
end

del(1)=B(1,1)/A(2,1);
for i=3:r+1
    for j=1:(size(A,2)-1)
        B(i,j)=B(i-2,j+1)-del(i-2)*A(i-1,j+1);
    end
    del(i-1)=B(i-1,1)/A(i,1);
end

P_1=0;
P_2=del(1);

Q_1=1;
Q_2=[1 gam(1)];

if r==1
	P=P_2;
	Q=Q_2;
end

for i=3:(r+1)
    if i>3
        P=del(i-1)*[1 zeros(1,(i-2))]+conv([1 0 0],P_1)+[0 gam(i-1)*P_2];
    else
        P=del(i-1)*[1 zeros(1,(i-2))]+[0 gam(i-1)*P_2];
    end
    Q=conv([1 0 0],Q_1)+[0 gam(i-1)*Q_2];
    P_1=P_2;P_2=P;
    Q_1=Q_2;Q_2=Q;
end
[P,Q]=tfdata(tf(P,Q)+D_fact,'v');
P=P.*(abs(P)>1e-6);
Q=Q.*(abs(Q)>1e-6);
RouthApprox=tf(P,Q)+D_fact;

