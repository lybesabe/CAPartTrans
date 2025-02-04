% Hopping Mechanism with Backscattering and Tunneling
clear; clc;
N = 10;
pd = 0.5;                   % probability to go down                                                             
pdi = 0.25;                 % probability to move diagonally (downwards)
pb = 0.5;                   % backscattering probability          
pt = 0.5;                   % tunneling probability                                                                     
iter = 30000;               % number of iterations per open-site vacancy

G = zeros(N, N);            % Grid as a matrix
Rx = 0:1:100;
Ry = zeros(1,101);

tic
for n=0:100 % open-site vacancy
    O=(N*N)-(N*N*(n/100));
    O=floor(O);

    a = 0;
    for k = 1:iter
        % construct the lattice (include the closed sites)
        m = 0;
        for i = 1:N
            for j = 1:N
                if m < O
                    G(i,j) = 1;
                    m = m+1;
                else
                    G(i,j) = 0;
                end
            end
        end
        
        % randomize the closed sites
        for i = 1:N
            for j = 1:N
                s = ceil(rand()*N);
                l = ceil(rand()*N);
                temp = G(i,j);
                G(i,j) = G(s,l);
                G(s,l) = temp;
            end
        end
        
        % place a particle on top of the lattice
        for j = 1:N
            w = ceil(rand()*N);
            if G(1,w) == 0
                G(1,w) = 2;
                break
            else
                j = 1;
            end
        end


        % Make particle travel
        j = w; i = 1;
        while i < N
            p = rand();
            q = rand();
            u = rand();
            x = ceil(rand()*100);
            t = (-1)^x;
            r = mod(j+t-1, N) + 1;
            c = mod(r+t-1, N) + 1;
            d = mod(j-t-1, N) + 1;
            if (p>=(1-pd) & G(i,j)==2 & G(i+1,j)==0)
                G(i+1,j) = 2;
                G(i,j) = 0;
                i = i + 1;
            elseif p>=pdi & p<1-pd & (G(i,j)==2 & G(i+1,r)==0)
                G(i+1, r) = 2;
                G(i,j) = 0;
                j = r;
                i = i + 1;
            elseif G(i,j) == 2 & G(i,r) == 0
                G(i,r) = 2;
                G(i,j) = 0;
                j = r;
            elseif q > 1-pt
                if p >= (1-pd) & i < (N-1) & G(i,j) == 2 & G(i+1,j) == 1 & G(i+2,j) == 0
                    G(i+2,j) = 2;
                    G(i,j) = 0;
                    i = i + 2;
                elseif p >= pdi & p < (1-pd) & i < (N-1) & G(i,j) == 2 & G(i+1,r) == 1 & G(i+2,c) == 0
                    G(i+2,c) = 2;
                    G(i,j) = 0;
                    i = i + 2;
                    j = c;
                elseif G(i,j) == 2 & G(i,r) == 1 & G(i,c) == 0
                    G(i,c) = 2;
                    G(i,j) = 0;
                    j = c;
                end
            elseif u > 1-pb
                if p >= (1-pd) & i > 1 & G(i,j) == 2 & G(i+1,j) == 1 & G(i-1,j) == 0
                    G(i-1,j) = 2;
                    G(i,j) = 0;
                    i = i - 1;
                elseif p >= pdi & p < (1-pd) & i > 1 & G(i,j) == 2 & G(i+1,r) == 1 & G(i-1,f) == 0
                    G(i+2,f) = 2;
                    G(i,j) = 0;
                    i = i - 1;
                    j = f;
                elseif G(i,j) == 2 & G(i,r) == 1 & G(i,f) == 0
                    G(i,f) = 2;
                    G(i,j) = 0;
                    j = f;
                end    
            else
                i = N;
            end
        end

        % Check if particle reached the bottom of the lattice
        if any(G(N, :) == 2)
            a = a + 1;
        end
    end % end of iteration

    q = a/iter; % percolation probability
    Ry(1,n+1) = q;
end % end of open site vacancy

toc
writematrix(Ry, 'perc_data_HM_10.csv')
plot(Rx,Ry);