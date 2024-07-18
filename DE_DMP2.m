function [ gbest, fitbest, iter] = DE_DMP2(q_min,q_max,q0i,E,e,accurat,ID_Run)
% output:
% gbest: bo nghiem tot nhat(rad)
% fitbest: gia tri loss cua ham muc tieu
% iter: so lan lap can thiet
% input:
% q_min,q_max: gioi han goc min,max (rad)
% q0i: gia tri q taij diem truoc do (rad)
% E: toa do diem cuoi (mm)
% e: diem thu e cua quy dao
% accurat: do chinh xac mong muon
% ID_Run: Chon thuat toan toi uu
D=7; %D: Number of variables;
iter_max=200; %so lan lap toi da
NP = 50;%quy mo dan do(so luong nghiem moi lan lap)
F_weight=0.5; Cr=0.9; I_bnd_constr=1;
FES_max= iter_max*NP;
%FES_max= 10000*D;
I_strategy=8;
%                1 --> DE/rand/1:                
%                2 --> DE/best/1:
%                3 --> DE/rand/2
%                4 --> DE/best/2
%                5 --> DE/Current to best/1
%                6 --> DE/Current to best/2 
%                7 --> DE/rand to best/1
%                8 --> DE/Current to rand/1  
% -----Initialize population and some arrays-------------------------------
%random varibales between the min and max values of the parameters
if ID_Run==1 || ID_Run==2 || ID_Run==3
    q_maxs=q_max;
    q_mins=q_min;
elseif ID_Run==4 || ID_Run==5 || ID_Run==6
    if e==2 %mo giong mien tim kiem o diem thu 2 cua quy dao
    q_max=q_max+pi/10;
    q_min=q_min-pi/10 ;      
    q_maxs=q0i+(q_max-q_min)/(50);%gioi han mien tim kiem quanh q0i tuy chinh theo tung quy dao
    q_mins=q0i-(q_max-q_min)/(50);
    else
    q_maxs=q0i+(q_max-q_min)/(10);
    q_mins=q0i-(q_max-q_min)/(10);

    
    end
end

          a=0;
 %    

    X=repmat(q_mins,NP,1)+rand(NP,D).*repmat(q_maxs-q_mins,NP,1); 
    
%      end
funevals = 0;      % number of function evaluations
f=ones(NP,1);

for i=1:NP
 f(i)=testfunctionDMP_arms(X(i,:),E,q0i,a);
end
funevals  = funevals + NP;   
%Find the best fitness in the population
[fitbest,index_min] = min(f);
gbest_iter = X(index_min,:); % best member of current iteration 
gbest = gbest_iter;            % best member ever
%------DE-Minimization---------------------------------------------
X_best    = zeros(NP,D);   % initialize gbestber  matrix
FVr_rot  = 0:1:NP-1; % rotating index array (size NP)
%create the best so far chart and average fitnesses chart.
Best_fitess=[];   Mean_fitness=[];
iter = 1;
%% Main loop of DE
while ((iter <= iter_max) && (fitbest >= accurat))
%while ((funevals <= FES_max) && (fitbest >= accurat))    
%while (iter <= iter_max)
%while (funevals <= FES_max)
%while (fitbest >= accurat)
    %ranking objective function    
    [f,f_index]=sort(f);
    X=X(f_index,:);
    %Choose random particles
    FVr_ind = randperm(4);    % index pointer array
    r1  = randperm(NP);       % shuffle locations of vectors
    FVr_rt  = rem(FVr_rot+FVr_ind(1),NP);   % rotate indices by ind(1) positions
    r2  = r1(FVr_rt+1);       % rotate vector locations
    FVr_rt  = rem(FVr_rot+FVr_ind(2),NP);
    r3  = r2(FVr_rt+1);       % rotate vector locations 
    FVr_rt  = rem(FVr_rot+FVr_ind(3),NP);
    r4  = r3(FVr_rt+1);       % rotate vector locations  
    FVr_rt  = rem(FVr_rot+FVr_ind(4),NP);
    r5  = r4(FVr_rt+1);       % rotate vector locations           
    X_r1 = X(r1,:);         % shuffled population 1
    X_r2 = X(r2,:);         % shuffled population 2
    X_r3 = X(r3,:);         % shuffled population 3
    X_r4 = X(r4,:);         % shuffled population 4
    X_r5 = X(r5,:);         % shuffled population 5
    X_best=repmat(gbest_iter,NP,1); % population filled with the best member
    %*-*-Mutation process-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      I_strategy=randi(8);
%  I_strategy=2;
    if (I_strategy == 1)  %1 --> DE/rand/1
        V = X_r1 + F_weight*(X_r2 - X_r3);   
        %FM_origin = X_r1;
    elseif (I_strategy == 2) %2 --> DE/best/1
        V = X_best + F_weight*(X_r1 - X_r2);
        %FM_origin = X_best; 
    elseif (I_strategy == 3) %3 --> DE/rand/2
        V = X_r1 + F_weight*(X_r2 - X_r3)+F_weight*(X_r4 - X_r5);               
        %FM_origin = X_r1;
    elseif (I_strategy == 4) %4 --> DE/best/2
        V = X_best + F_weight*(X_r1 - X_r2)+F_weight*(X_r3 - X_r4);
        %FM_origin = X_best;
    elseif (I_strategy == 5) %5 --> DE/Current to best/1
        V = X + F_weight*(X_best - X)+F_weight*(X_r1 - X_r2);
        %FM_origin = X;
    elseif (I_strategy == 6) %6 --> DE/Current to best/2 
        V = X + F_weight*(X_best - X)+F_weight*(X_r1 - X_r2)+F_weight*(X_r3 - X_r4);
        %FM_origin = X;  
    elseif (I_strategy == 7) %7 --> DE/rand to best/1
        V = X_r1 + F_weight*(X_best - X_r1)+F_weight*(X_r2 - X_r3);
        %FM_origin = X_r1;
    elseif (I_strategy == 8) %8 --> DE/Current to rand/1
        V = X + F_weight*(X_r1 - X)+F_weight*(X_r2 - X_r3);
        %FM_origin = X;
    end
    %*-*-*Crossover Process*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    D_index=repmat(1:D,NP,1);
    Jrand=randi(D,NP,1); % Interger random number matrix[NPx1]  in [1, D]   
    Jrand=repmat(Jrand,1,D); %copy collom 1 to the others
    FM_mui = (rand(NP,D) < Cr)|(Jrand==D_index);  % all random numbers < F_CR are 1, 0 otherwise
    FM_mpo = FM_mui < 0.5;    % inverse mask to FM_mui 
    U = X.*FM_mpo + V.*FM_mui  ;    
    %=====Only use this if boundary constraints are needed==================
    if (I_bnd_constr == 1)
%         gioi han nghiem tim kiem trong mien tim kiem
        UlowerIndex=U<repmat(q_mins,NP,1); UupperIndex=U>repmat(q_maxs,NP,1);
        U=repmat(q_mins,NP,1).*UlowerIndex+U.*(~UlowerIndex);    
        U=repmat(q_maxs,NP,1).*UupperIndex+U.*(~UupperIndex);
%         gioi han tim kiem theo gia tri goc max,min
%         UlowerIndex=U<repmat(q_min,NP,1); UupperIndex=U>repmat(q_max,NP,1);
%         U=repmat(q_min,NP,1).*UlowerIndex+U.*(~UlowerIndex);    
%         U=repmat(q_max,NP,1).*UupperIndex+U.*(~UupperIndex);
        
    end
    %=====End boundary constraints==========================================
    %-----Optional parent+child selection-----------------------------------------
    %-----Select which vectors are allowed to enter the new population------------
    for k=1:NP
        Score = testfunctionDMP_arms(U(k,:),E,q0i,a);   % check cost of competitor
        funevals  = funevals + 1;
        if (Score<=f(k))   
            X(k,:) = U(k,:);                    % replace old vector with new one (for new iteration)
            f(k)   = Score;                      % save value in "cost array"
            %----we update GlobalMin only in case of success to save time-----------
            if (Score<=fitbest)   
                fitbest = Score;                    % new best value
                gbest = U(k,:);                 % new best parameter vector ever
            end
        end      
    end % end of "for k = 1:NP"
    gbest_iter = gbest;       % freeze the best member of this iteration for the coming 
    %----Output section----------------------------------------------------------    
    % formatSpec ='DE: iter: %d and fitbest: %0.5e /n';
    % fprintf(formatSpec,iter,fitbest);
    % for n=1:D
    %    fprintf(1,'gbest(%3d) = %-0.15f\n',n,gbest(n));
    % end
    Best_fitess=[Best_fitess fitbest];
    Mean_fitness=[Mean_fitness mean(f)];
    iter = iter + 1;
end 
end

