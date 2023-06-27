
%% Problem Definition
function F2 = mainccring_2()

fobj=@(x) MOP2(x);      % Cost Function

nVar=2;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix
VarMin=[-4,-4];          % Lower Bound of Variables
VarMax=[4,4];          % Upper Bound of Variables
VarMin1=[-4,-4];          % Lower Bound of Variables
VarMax1=[4,4];          % Upper Bound of Variables
% VarMin= [0,0,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2];
% VarMax = [1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2];
% VarMin1= [0,0,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2];
% VarMax1 = [1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2];
% VarMin=[0,0];          % Lower Bound of Variables
% VarMax=[1,1];          % Upper Bound of Variables

%%NCCLA parameters
MaxIt=500;      % Maximum Number of Iterations
RPprob=0.9; %reinforcement learning
SLprob=0.99; %social learning
VSLprob=0.99; %vector or horizontal333
P1prob=0.95; %first parent
Taeprob=0.3; % trial and error
ifmin=0.0005; %
ifmax=0.02;
iter=1;
npop=100;% Population Size
sswarm=1:1:npop;
rep_size=100;
%lb=-5;
%ub=10;


d=nVar;
%funnum=14;
% BestSol.Cost=inf;



%% Initialization

disp('Initialization ...');

empty_individual.position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance_obj=[];
empty_individual.CrowdingDistance_dec=[];
empty_individual.CrowdingDistance=[];

 empty_individual.pbest=[];
 empty_individual.nbest=[];
 empty_individual.PBA=[];
 empty_individual.NBA=[];


pop=repmat(empty_individual,npop,1);
rep=repmat(empty_individual,rep_size,1);
%1. start build first population

for t=1:npop
pop(t).position=initialization(nVar,VarMax,VarMin);
    % Evaluate the fitness of the new solution
%pop(t).Cost=CostFunction(pop(t).position); 
pop(t).Cost=fobj(pop(t).position); 

end

  [pop F]=NonDominatedSorting(pop);

% Calculate Crowding Distance
 pop=CalcCrowdingDistance(pop,F);

% Sort Population
 [pop F]=SortPopulation(pop);
for t=1:npop
     pop(t).PBA=pop(t);
end
%dtermine best neighbor
for i=1:npop
    if i==1
%         pop(i,1).NBA=NEdetermine(pop(npop).PBA,pop(randi(npop),1).PBA,pop(i+1,1).PBA);
                pop(i,1).NBA=NEdetermine(pop(npop).PBA,pop(1,1).PBA,pop(i+1,1).PBA);

    end
    
    if i==npop
%          pop(i,1).NBA=NEdetermine(pop(npop-1,1).PBA,pop(randi(npop),1).PBA,pop(1,1).PBA);
                  pop(i,1).NBA=NEdetermine(pop(npop-1,1).PBA,pop(npop,1).PBA,pop(1,1).PBA);


    end
    
    if i>1 & i<npop
%          pop(i,1).NBA=NEdetermine(pop(i-1,1).PBA,pop(randi(npop),1).PBA,pop(i+1,1).PBA);
                  pop(i,1).NBA=NEdetermine(pop(i-1,1).PBA,pop(i,1).PBA,pop(i+1,1).PBA);


    end
end 
Ffirst=pop(F{1});
ranks=[pop.Rank];
idmaxrank=find(ranks==max(ranks));
if idmaxrank>1
      idxworst=randi([1 size(idmaxrank,1)]);
      WorstSol.Cost=pop(idmaxrank(1,idxworst),1).Cost;
      WorstSol.position=pop(idmaxrank(1,idxworst),1).position;
else
    WorstSol.Cost=pop(idmaxrank(1,1),1).Cost;
    WorstSol.position=pop(idmaxrank(1,1),1).position;
end
%2.determine parents//juveniles

%idxbest=randi([1 size(Ffirst,1)],2,1);
idxbest=randperm(size(Ffirst,1));

BestSol=Ffirst(idxbest(1,1),1);
% BestSol.position=Ffirst(idxbest(1,1),1).position;
%3.main loop//social,asocial,learning phase
rep(1:size(Ffirst,1),:)=Ffirst;
%2.determine parents//juveniles

%3.main loop//social,asocial,learning phase
 while iter<=MaxIt
 a(iter)=0.5*(log((1+(iter/MaxIt))/(1-(iter/MaxIt))));
 vb= (a(iter)+a(iter)).*rand(1,1)-a(iter);
 a1=2-(2-0).*(((exp(iter/MaxIt)-1)/(exp(1)-1))^0.2);
 A2=2.*a1.*rand-a1;
 juv=pop;

for i=1:size(juv,1)
    a7=randperm(3);
    p1=pop(i,1).NBA(a7(1,1));
    p2=pop(i,1).NBA(a7(1,2));
    ccc=[juv.Cost];
    A3=tanh(((mean(p1.Cost)+mean(p2.Cost))/2)-(mean(ccc(:,1))+mean(ccc(:,2)))/2);
    a2=size(pop(i,1).PBA,1);
    a2=randperm(a2);
    pbest=pop(i,1).PBA(a2(1,1));

    for j=1:d
        rand1=rand;
         if rand1<=SLprob
%         if rand1<=A3

            rand2=rand;
            if rand2<=VSLprob
                rand4=rand;
                if rand4<=P1prob
                    pop(i).position(j)=p1.position(j);   
                else
                    pop(i).position(j)=p2.position(j);
                end
            else
                r = randperm(npop);
                if r(1,1)~=i
                pop(i).position(j)=pop(r(1,1)).position(j); 
                else
                pop(i).position(j)=pop(r(1,2)).position(j); 

                end
            end          
         else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             yk=rand(1,nVar);
%              for z1=1:nVar
%             
%              ykn(1,z1)=4.*yk(1,z1).*(VarMax1(1,z1)-VarMin1(1,z1));
%              
%              end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rand3=rand;
            if rand3<=0.75
                 pop(i).position(j)=rand(1,1).*(VarMax(1,j)-VarMin(1,j))+VarMin(1,j);
%                 juv(i).position(j)=ykn(j).*(VarMax1(1,j)-VarMin1(1,j))+VarMin1(1,j);

            else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%                 cau=cauchyrnd(0, 1, nVar);
%                 cau=cau(1,:);
%                 cau=rescale(cau,-1,1);
%                 juv(i).position=juv(i).position+cau;
                 
%                  pop(i).position(j)=cau(1,j).*pbest.position(j);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  pop(i).position(j)=pbest.position(j);

            end
            
        end
    end
    Flag4ub=juv(i).position>VarMax;
    Flag4lb=juv(i).position<VarMin;
    juv(i).position=(juv(i).position.*(~(Flag4ub+Flag4lb)))+VarMax.*Flag4ub+VarMin.*Flag4lb;
    iff=ifmin+((ifmax-ifmin)/MaxIt)*iter;
    pop(i).position=RLJ(pop(i).position,juv(i).position,juv,i,RPprob,iff);
    Flag4ub=juv(i).position>VarMax;
    Flag4lb=juv(i).position<VarMin;
    juv(i).position=(juv(i).position.*(~(Flag4ub+Flag4lb)))+VarMax.*Flag4ub+VarMin.*Flag4lb;
    pop(i).Cost=fobj(pop(i).position); 
   [pop(i,1).NBA(a7(1,1)).position,pop(i,1).NBA(a7(1,2)).position]=RLP(p1.position,p2.position,juv);

end

% merge
pop=[pop  
    juv];%%merge

%Non-dominated sorting
[pop F3]=NonDominatedSorting(pop);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F3);

% Sort Population
[pop F3]=SortPopulation(pop);
pop=pop(1:npop);
%%update PBA
for i=1:npop
    tedadpbest=size(pop(i,1).PBA,1);
    
pop(i,1).PBA(tedadpbest+1,1)=pop(i,1);

 [pop(i,1).PBA F6]=NonDominatedSorting(pop(i,1).PBA);

% Calculate Crowding Distance
pop(i,1).PBA=CalcCrowdingDistance(pop(i,1).PBA,F6);

% Sort Population
[pop(i,1).PBA F6]=SortPopulation(pop(i,1).PBA);
if size(pop(i,1).PBA,1)>5
 pop(i,1).PBA=pop(i,1).PBA(1:5,1);
end
end

%%Update NBA
for i=1:npop
    if i==1
         pop(i,1).NBA=NEdetermine(pop(npop,1),pop(1,1),pop(2,1));
%                 pop(i,1).NBA=NEdetermine(pop(npop,1).PBA,pop(1,1).PBA,pop(i+1,1).PBA);

    end
    
    if i==npop
         pop(i,1).NBA=NEdetermine(pop(npop-1,1),pop(npop,1),pop(1,1));
%                   pop(i,1).NBA=NEdetermine(pop(npop-1,1).PBA,pop(npop,1).PBA,pop(1,1).PBA);


    end
    
    if i>1 & i<npop
          pop(i,1).NBA=NEdetermine(pop(i-1,1),pop(i,1),pop(i+1,1));
%                   pop(i,1).NBA=NEdetermine(pop(i-1,1).PBA,pop(i,1).PBA,pop(i+1,1).PBA);


    end
    pop(i,1).NBA=pop(i,1).NBA(1:3,1);
end


 %Non-dominated sorting
[pop F]=NonDominatedSorting(pop);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F);

% Sort Population
[pop F]=SortPopulation(pop);


for i=1:npop
    %%update NBA
 [pop(i,1).NBA F5]=NonDominatedSorting(pop(i,1).NBA);

% Calculate Crowding Distance
pop(i,1).NBA=CalcCrowdingDistance(pop(i,1).NBA,F5);

% Sort Population
[pop(i,1).NBA F5]=SortPopulation(pop(i,1).NBA);
pop(i,1).NBA=pop(i,1).NBA(1:3,1);

end
p=0;
for i=1:npop
    r=size(pop(i,1).NBA,1);
    if i==1
     r=r+p;
z3(i:r,1)=pop(i,1).NBA(1:end,1);
    else
     r=r+p-1;

z3(p:r,1)=pop(i,1).NBA(1:end,1);
    end
p=r+1;

end

[z3 F1]=NonDominatedSorting(z3);

% Calculate Crowding Distance
z3=CalcCrowdingDistance(z3,F1);

% Sort Population
[z3 F1]=SortPopulation(z3);

% [pop F]=NonDominatedSorting(pop);
% % 
% % % Calculate Crowding Distance
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  pop=CalcCrowdingDistance(pop,F);
% % % % % 
% % % % % % Sort Population
% %  [pop F]=SortPopulation(pop);

% % % % Store F1
 if iter>1
     F1=[z3(F1{1});pop(F{1})];
    
    [F1 FN]=NonDominatedSorting(F1);

% Calculate Crowding Distance
F1=CalcCrowdingDistance(F1,FN);

% Sort Population
[F1 FN]=SortPopulation(F1);
F1=F1(FN{1});
if size(F1,1)>100
F1=F1(1:100,1);
end
 else
     F1=pop(F{1});
 end
 Fupdate=pop(F{1});
    
idxbest=randperm(size(Fupdate,1));
g=[pop.Rank];
if size(find(g)==1,2)==1
Fsecond=pop(F{2});
secid=randi(size(Fsecond,1));
idxbest(1,2)=Fsecond(secid);
end

%BestSol=Ffirst(idxbest(1,1),1);
BestSol=Fupdate(idxbest(1,1),1);

ranks=[pop.Rank];
idmaxrank=find(ranks==max(ranks));
if idmaxrank>1
      idxworst=randi([1 size(idmaxrank,1)]);
      WorstSol.Cost=pop(idmaxrank(1,idxworst),1).Cost;
      WorstSol.position=pop(idmaxrank(1,idxworst),1).position;
else
    WorstSol.Cost=pop(idmaxrank(1,1),1).Cost;
    WorstSol.position=pop(idmaxrank(1,1),1).position;
end
Fupdate=[Fupdate
         Ffirst];

     
       % Non-Dominated Sorting
    [Fupdate Ffinal]=NonDominatedSorting(Fupdate);

    % Calculate Crowding Distance
    Fupdate=CalcCrowdingDistance(Fupdate,Ffinal);

    % Sort Population
    [Fupdate Ffinal]=SortPopulation(Fupdate);
     if size(Fupdate,1)>100
         Fupdate=Fupdate(1:rep_size);
                % Non-Dominated Sorting
    [Fupdate Ffinal]=NonDominatedSorting(Fupdate);

    % Calculate Crowding Distance
    Fupdate=CalcCrowdingDistance(Fupdate,Ffinal);

    % Sort Population
    [Fupdate Ffinal]=SortPopulation(Fupdate);

        F2=Fupdate(Ffinal{1});
     else
        F2=Fupdate(Ffinal{1});
     end
%       F=[];
 disp(['Iteration ' num2str(iter) ': Number of F1 Members = ' num2str(numel(F2))]) 
%  end
 
% % %  display(['At iteration ', num2str(iter), ' the elite fitness is ', num2str(BestSol.fitness)])
%  gbest(1,iter)=BestSol.fitness;
% Plot F1 Costs
%    if iter>1
%         figure(1);
%     PlotCosts(FNN);
%    else 


    PlotCosts(F2);

% A = imread('my_img.png');
% for ii = 1:N   
%    imwrite(A,strcat('my_new',num2str(ii),'.tiff'));
% end
%     exportgraphics(f,'barchart.png','Resolution',300)
% print(gcf(iter), '-dtiff', 'myfigure.tiff');
if iter==1||iter==10||iter==25||iter==40||iter==70||iter==100||iter==500
print(gcf, '-dtiff', strcat('my_new_MO_ncc_kham',num2str(iter),'.tiff'));
end

%     end
    
%     Fss{1}=[Fss{1}; F1];
iter=iter+1;
 end
 
  quite
end
 