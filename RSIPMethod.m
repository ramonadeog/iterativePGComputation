function [S,err] = RSIPMethod(B,Ns,Nr,C,numIter,txInd,Txx)
%RSIPMETHOD 
%S: State Matrix at convergence
%B: scatterer to scatterer matrix for the entire graph
%Ns: number of scatterer per room
%Nr: number of rooms
%C: Nr by 1 cell array containing list of neighbours for each room

Aii = zeros(Nr,Ns,Ns);
for ii = 1:Nr
    Aii(ii,:,:)=(eye(Ns)-B((ii-1)*Ns+1:ii*Ns,(ii-1)*Ns+1:ii*Ns))\eye(Ns);
end

SS = zeros(Nr,Ns,numIter+1);
for ii = 1:numIter
    for jj = 1:Nr
        neighB = cell2mat(C(jj));
        if jj==txInd
            Svalue = Txx;
        else
            Svalue = 0;          
        end
        for uu = 1:length(neighB)
            Svalue = Svalue+B((jj-1)*Ns+1:jj*Ns,(neighB(uu)-1)*Ns+1:neighB(uu)*Ns)*squeeze(SS(neighB(uu),:,ii)).';
        end
        SS(jj,:,ii+1) = squeeze(Aii(jj,:,:))*Svalue;
    end
    err(ii) = norm(squeeze(SS(:,:,ii+1)-SS(:,:,ii)))/norm(squeeze(SS(:,:,ii)));
end
S=squeeze(SS(:,:,numIter+1));



end

