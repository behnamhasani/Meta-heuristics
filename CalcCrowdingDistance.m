function pop=CalcCrowdingDistance(pop,F)

    nF=numel(F);
     n_var=size(pop(1).position,2);
% 
%     for i = 1 : n_var
%        crowd_dist_var(:,1) = crowd_dist_var(:,1) + y(:,n_obj + n_var + 1 + i);
%     end
%     
%     crowd_dist_var = [];
%         crowd_dist_var(:,1) = zeros(length(F(front).f),1);
        
    for k=1:nF
        
        
        Costs=[pop(F{k}).Cost];
        
        nObj=size(Costs,1);
        
        n=numel(F{k});
% %         tt=zeros(1,n);
% %         tt=F{k};
% %         
% %         for i=1:n
% %            for j = 1 : n_var
% %                 pos(i,:) = pop(tt(1,i)).position;
% %            end
% %         end
        
        d=zeros(n,nObj);
        
        for j=1:nObj
            
            [cj so]=sort(Costs(j,:));
            
            d(so(1),j)=inf;
            
            for i=2:n-1
                
                d(so(i),j)=abs(cj(i+1)-cj(i-1))/abs(cj(1)-cj(end));
                
            end
            
            d(so(end),j)=inf;
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%decision space
%         d1=zeros(n,n_var);
%           for j=1:n_var
%             
%             [cj1 so1]=sort(pos(:,j));
%             
%             d1(so1(1),j)=inf;
%             
%             for i=2:n-1
%                 
%                 d1(so1(i),j)=abs(cj1(i+1)-cj1(i-1))/abs(cj1(1)-cj1(end));
%                 
%             end
%             
%             d1(so1(end),j)=inf;
%             
%           end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        avg_distdec=0;
        avg_distobj=0;
%         for i=1:n
%             
%              pop(F{k}(i)).CrowdingDistance_dec=sum(d1(i,:));
%                t2=d1(i,:);
% %             pop(F{k}(i)).CrowdingDistance_dec=sum(t2(isfinite(t2)));
%             avg_distdec=avg_distdec+sum(t2(isfinite(t2)));
%             
% %         end
%         
% %         for i=1:n
%             
%             pop(F{k}(i)).CrowdingDistance_obj=sum(d(i,:));
%             t3=d(i,:);
%             avg_distobj=avg_distobj+sum(t3(isfinite(t3)));
%             
%         end
%          avg_distdec=avg_distdec/n;
%          avg_distobj=avg_distobj/n;
         for i=1:n
%             if (pop(F{k}(i)).CrowdingDistance_dec>avg_distdec)||(pop(F{k}(i)).CrowdingDistance_obj>avg_distobj)
%                pop(F{k}(i)).CrowdingDistance=max(pop(F{k}(i)).CrowdingDistance_dec,pop(F{k}(i)).CrowdingDistance_obj);
%             else
%                pop(F{k}(i)).CrowdingDistance=min(pop(F{k}(i)).CrowdingDistance_dec,pop(F{k}(i)).CrowdingDistance_obj);
%             end

          pop(F{k}(i)).CrowdingDistance=sum(d(i,:));  
        end
         
        
    end


end

