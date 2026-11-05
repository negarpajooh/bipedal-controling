function f=mpcsetup(A,B,C,D,p,m,Q,R,x0,u0)


% Input and output check
error(nargchk(8,10,nargin));
error(nargoutchk(1,1,nargout));

% dimension of the system
nx=size(A,1);
[ny,nu]=size(D);

% precalculate the gain matrices of the online controller
[K,Pr] = predmat(A,B,C,D,p,m,Q,R);
% only the first instance is required
K1 = K(1:nu,:);
Pr = Pr(1:nu,:);

% default initial conditions
if nargin<9
    x0=zeros(nx,1);
end
if nargin<10
    u0=zeros(nu,1);
end

% the online controller
f = @ssmpc;
    function u=ssmpc(y,r)
        % update state
        x0 = A*x0 + B*u0;
        % get the reference
        nr=size(r,2);
        if nr>=p
            ref=r(:,1:p);
        else
            ref=[r r(:,end+zeros(p-nr,1))];
        end
        % model-plant mismatch
        offset = (y(:) - C*x0);
        % update optimal control
%         u0 = -K1*x0 + Pr*reshape(bsxfun(@minus,ref,offset),[],1);
        u0 = -K1*x0 + Pr*reshape(ref-offset(:,ones(p,1)),[],1);
        % controller output
        u = u0;
    end
end


function [K,Pr] = predmat(A,B,C,D,hP,hM,Q,R)
 
%%%% Initialise
[nx,nu]=size(B);
ny = size(C,1);

Gss = inv([A-eye(nx) B;C D]);
M1 = Gss(nx+1:end,nx+1:end);

Px=zeros(hP*ny,nx);  Pu=zeros(hP*ny,hP*nu);  P=C;
% L=eye(ny*hP);

%%%% Use recursion to find predictions
for i=1:hP;
   Puterm = P*B;
   for j=i:hP;
         vrow=(j-1)*ny+1:j*ny;
         vcol=(j-i)*nu+1:(j-i+1)*nu;
         Pu(vrow,vcol)=Puterm;
   end
   P=P*A;
   vrow=(i-1)*ny+1:i*ny;
   Px(vrow,1:nx) = P;
end
P=Px;
M=zeros(nu*hM,ny*hP);
for i=1:hM,
    I=(i-1)*nu+1:i*nu;
    J=(i-1)*ny+1:i*ny;
    M(I,J)=M1;
end

%%% Sum last columns of H
v=(hM-1)*nu;
HH=zeros(hP*ny,nu);
for k=1:nu;
    HH(:,k) = sum(Pu(:,v+k:nu:end),2);
end
H = [Pu(:,1:(hM-1)*nu),HH];

S = H'*Q*H+R;
S = (S+S')/2;
X1 = H'*Q*P;
X2 = -H'*Q-R*M;

%%%% Unconstrained control law
Mi=inv(S);
K = Mi*X1;
Pr = -Mi*X2;
end