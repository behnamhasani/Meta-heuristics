function z=MOP2(x)

        n=numel(x);
    sz1=0;
    sz2=0;
    for i=1:n
        sz1=sz1+(x(1,i)-1/sqrt(n))^2;
        sz2=sz2+(x(1,i)+1/sqrt(n))^2;

%     z1=1-exp(-sum((x-1/sqrt(n)).^2));
%     
%     z2=1-exp(-sum((x+1/sqrt(n)).^2));
    end
         z1=1-exp(-sz1);
%     
     z2=1-exp(-sz2);
    z=[z1 z2]';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      y = [0; 0; 0];
%     x=x.';
%     [dim, num]  = size(x);
%     Y           = zeros(dim,num);
%     Y(3:dim,:)  = (x(3:dim,:) - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
%     tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
%     tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
%     tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
%     y(1)      = cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2);
%     y(2)      = cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + 2.0*tmp2/size(5:3:dim,2);
%     y(3)      = sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end