function PlotCosts(pop)

    Costs=[pop.Cost];
    figure(1)
    plot3(Costs(1,:),Costs(2,:),Costs(3,:),'ro','MarkerSize',8,'MarkerFaceColor','r');
    xlabel('1st Objective');
    ylabel('2nd Objective');
    grid on;
    
    for i=1:size(pop,1)
          for j=1:3
            ngpos(i,j)=pop(i).position(1,j);
          end
    end
 figure(2)
    plot3(ngpos(:,1),ngpos(:,2),ngpos(:,3),'r*','MarkerSize',8,'MarkerFaceColor','r');
    xlabel('1st Objective');
    ylabel('2nd Objective');
    grid on;
    
end