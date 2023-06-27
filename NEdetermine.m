function z=NEdetermine(Nghabli,khodesh,Nbadi)
z=[Nghabli;
   khodesh;
   Nbadi];
[z F]=NonDominatedSorting(z);

% Calculate Crowding Distance
 z=CalcCrowdingDistance(z,F);

% Sort Population
 [z F]=SortPopulation(z);
    
end