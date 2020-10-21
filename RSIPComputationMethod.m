%Iterative State Matrix Computation for Multi-room Environment
%R. O. Adeogun, A. Bharti and T. Pedersen, "An Iterative Transfer Matrix 
%Computation Method for Propagation Graphs in Multi-Room Environments," 
%in IEEE Antennas and Wireless Propagation Letters.
%Feb. 2019
%======================Parameters========================================
clear all
clc
Nt = 1;                              %Number of Transmitters/antennas
Nr = 1;                              %Number of receivers/rx antennas
Nm = 10;                             %Number of rooms
Ns = 10;                             %Number of scatterers per room(Can be 
g = 0.52;                            %Reflection co-efficient
Pvis =0.80;                          %Probability of visibility

%=======Geometry Parameters=============================================
%roomOriginLoc = [0 0 0;0 4 0;0 8 0;0 12 0]; %Plan S1
%roomOriginLoc =[0 4 0; 3 4 0; 0 0 0; 3 0 0]; %Plan S2
%roomOriginLoc =[0 4 0;0 0 0;3 0 0;6 0 0;9 0 0;12 0 0; 15 0 0]; % Plan S3
roomOriginLoc = [0 4 0;0 0 0;3 4 0;3 0 0;6 4 0;6 0 0;9 4 0; 9 0 0;12 4 0....
                 ;12 0 0;15 4 0;15 0 0];
txInd = 1;     rxInd = 1;             %Index of transmit and receive rooms

%CNeigh = {[2 3];[1]};
%CNeigh = {[2];[1 3];[2 4];[3]};
%CNeigh = {[2 3 4 5 6 7];[1 3];[1 2 4];[1 3 5];[1 4 6];[1 5 7];[1 6]};
CNeigh = {[2 7 8];[1 3 6 8];[2 4 6];[3 5 6];[4 6 10];[2 3 4 5 8 9 10];...
               [1 8];[1 2 6 7 9];[6 8 10];[5 6 9]}; % Neighbouring rooms
              
Adj = zeros(Nm,Nm);
for ii = 1:Nm
    Adj(ii,cell2mat(CNeigh(ii)))=1;
end
Adj(logical(eye(size(Adj))))= 1;
AdjFull = kron(Adj,ones(Ns,Ns));
Rmm = [3 4 3; 2 2 3;3 2 3;3 2 3;6 5 3;6 2 3;3 4 3;2 6 3;3 4 3;3 4 3];

%==================================
c = 3e8;
fmin = 58e9; fmax = 62e9;
fc = 60e9;
Deltaf = 5e6;
f = fmin:Deltaf:fmax;
Nf = length(f);
Deltat = 1/(fmax-fmin); 
Taxis = (0:Nf-1)*Deltat*1e9;
N = Ns*Nm+Nr+Nt;
Nss = Nm*Ns;
wp = 0.2; pl=1; 
%numIter = 3:2:40;
numIter = 5; timePG = 0; timeISCM = 0;
PLoss = wp*ones(Nm,Nm);
PLoss(logical(eye(size(PLoss))))= 1;
PL = kron(PLoss,ones(Ns,Ns));
cc = 1; cfail = 0;
Count = 1; PP1=0; PP2 = 0; PP3 = 0;
Pt = 1; freqCorr = 0; freqCorrw=0; freqCorrwa=0; freqCorrii=0;
time1 = 0; time2=0;  errorr=0;
while cc <= Count
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%Adjacency Matrix Generation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E = zeros(N,N);
    E(Nt+Nr+1:N,Nt+Nr+1:N)=(arrayfun(@(z)sum(z >= cumsum([0 1-Pvis Pvis])), rand(Nm*Ns))-1).*AdjFull;
    E(logical(eye(size(E))))= 0;
    E(Nt+Nr+(txInd-1)*Ns+1:Nt+Nr+txInd*Ns,1:Nt) = (arrayfun(@(z)sum(z >=cumsum([0 1-Pvis Pvis])), rand(1,Ns))-1); 
    E(Nt+1:Nt+Nr,Nt+Nr+(rxInd-1)*Ns+1:Nt+Nr+rxInd*Ns) = (arrayfun(@(z)sum(z >=cumsum([0 1-Pvis Pvis])), rand(1,Ns))-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Scatterer Placement Within (or on room walls) Rooms 
    Rm = [];
    for ii = 1:Nm
        Rm = [Rm;scattererPlacement(Ns,Rmm(ii,:),2)+roomOriginLoc(ii,:)];
    end
    txLoc =scattererPlacement(Nt,Rmm(txInd,:),1)+roomOriginLoc(txInd,:);   %Transmitter in Room txInd
    
    rxLoc =scattererPlacement(Nr,Rmm(rxInd,:),1)+roomOriginLoc(rxInd,:);  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Rc = [txLoc; rxLoc; Rm];
    d = zeros(N,1);
    Psi_up = triu(0+2*pi*rand(N),1);
    Psi_e = diag(d)+Psi_up+Psi_up.';
    for ic = 1:Nf
        A = genWeightAdjacency(g,Nm*Ns,E,Nr,Nt,f(ic),Rc,Psi_e);
        B= A(Nt+Nr+1:Nt+Nr+Nss,Nt+Nr+1:Nt+Nr+Nss);
        B= PL.*B;
        if max(abs(eig(B)))>=1
            Hpg = 0;
            break
        end
        %%%%%%%%%%%%%%%%Original PG Computation=========================
        T = A(Nt+Nr+1:Nt+Nr+Nss,1:Nt);
        R = A(Nt+1:Nt+Nr,Nt+Nr+1:Nt+Nr+Nss); 
        Hpg(:,:,ic) = R*((eye(Nss)-B)\T); 
        %===============Recursive Inter-room Propagation Method==========
        Txx = T((txInd-1)*Ns+1:txInd*Ns,:);
        Rxx = R(:,(rxInd-1)*Ns+1:rxInd*Ns);
        S = RSIPMethod(B,Ns,Nm,CNeigh,numIter,txInd,Txx);
        Hiscm(:,:,ic) = Rxx*S(rxInd,:).';
        %========================================================
 
    end
    if Hpg == 0
        cfail = cfail+1;
    else
        for jj=1:Nr
            for kk=1:Nt
                hpg(jj,kk,:) =(abs(ifft(squeeze(Hpg(kk,jj,:))))); 
                hiscm(jj,kk,:) =(abs(ifft(squeeze(Hiscm(kk,jj,:)))));
            end
        end
        pp1= (squeeze(mean(mean(abs(hpg).^2,1),2)));
        pp2 = (squeeze(mean(mean(abs(hiscm).^2,1),2)));
        PP1 = pp1+PP1;
        PP2 = pp2+PP2;
        cc = cc+1
    end   
end
PP1 = PP1/Count; 
PP2 = PP2/Count;
figure(1);
plot(Taxis,10*log10(PP1),'linewidth',2); xlim([0 150])
hold on
plot(Taxis,10*log10(PP2),'--','linewidth',1); xlim([0 150])
xlabel('Delay [ns]'); ylabel('Power [dB]')