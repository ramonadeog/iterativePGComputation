function [A] = genWeightAdjacency(g,Ns,E,Nr,Nt,f,Rc,Psi_e)
%GENWEIGHTADJACENCY 
%Generate weighted adjacency matrix A(f) for a unipolarized graph 
%with Nr receiver(s) and Nt transmitter(s) at freq. f
%rng('default')

N = Ns+Nr+Nt;
%E(1:Nt,:) = 0;                                %Transmitter is a source
%E(:,Nt+1:Nt+Nr) = 0;                          %Receiver is a sink
%E(Nt+1:Nt+Nr,1:Nt) = 1;                       %Include Direct path
%E(logical(eye(size(E))))= 0;                  %No loops
%
% Rscx = Rgrid(1)*rand(Ns,1);
% Rscy = Rgrid(2)*rand(Ns,1);
% Rscz = Rgrid(3)*rand(Ns,1);
% Rsc =[Rscx(:) Rscy(:) Rscz(:)];

dis_e = pdist2(Rc,Rc);  %Distance on all possible edges
tau_e = dis_e/3e8;      %Delay on all possible edges


Dtau = tau_e(Nt+1:Nt+Nr,1:Nt);
Ttau = tau_e(Nt+Nr+1:Nt+Nr+Ns,1:Nt);
Tcard = (find(Ttau));
Rtau = tau_e(Nt+1:Nt+Nr,Nt+Nr+1:Nt+Nr+Ns);
Rcard = (find(Rtau));
Btau = tau_e(Nt+Nr+1:Nt+Nr+Ns,Nt+Nr+1:Nt+Nr+Ns);   
Bcard = (find(Btau));
deg_B = zeros(Ns,Ns);
Bg = zeros(Ns,Ns);
BEdge = E(Nt+Nr+1:Nt+Nr+Ns,Nt+Nr+1:Nt+Nr+Ns);
for jj  = 1:Ns
    deg_B(:,jj) = ones(1,Ns)*length(find(Btau(:,jj)));
end
 %mus = mean(Btau(Bcard))*1e9;
 %g = db2pow(rho*mus/2);
Bg(Bcard) = (g^2.*BEdge(Bcard))./deg_B(Bcard);
Ag = zeros(N,N);
%Compute Edge gain matrix.
Dg = 1./(4*pi*f.*Dtau);
Ag(Nt+1:Nt+Nr,1:Nt) = Dg;
Tg = zeros(Ns,Nt);
Tg(Tcard) = (2*Nt*length(Tcard)./(4*pi*f*sum(sum(Ttau)))).*(1./(Ttau(Tcard).^(2).*sum(sum(Ttau(Tcard).^(-2)))));
Ag(Nt+Nr+1:Nt+Nr+Ns,1:Nt) = sqrt(Tg);
Rg = zeros(Nr,Ns);
Rg(Rcard) = (2*Nr*length(Rcard)./(4*pi*f*sum(sum(Rtau)))).*(1./(Rtau(Rcard).^(2).*sum(sum(Rtau(Rcard).^(-2)))));
Ag(Nt+1:Nt+Nr,Nt+Nr+1:Nt+Nr+Ns) = sqrt(Rg);
Ag(Nt+Nr+1:Nt+Nr+Ns,Nt+Nr+1:Nt+Nr+Ns) = sqrt(Bg);
A = E.*Ag.*exp(1j*(Psi_e-2*pi*tau_e*f));

end

