x=comp(:,2);
y=comp(:,1);

for i=1:length(x)
    
        if(y(i)==y(i+1))
        
            x(i+1)= x(i)+x(i+1);
            x(i)=0;
            end
        end
     
        
            
