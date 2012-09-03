 
SortedData = sortrows(data,[2 1]);

 for i = 1:41   
        for k = 1:7
            for j = 1:101
                Liang(i,j,k)=SortedData((i-1)*101 + j,k);
                Liang(i,102,k)=SortedData((i-1)*101 + 1,k);                
            end
        end 
 end
 for i = 1:41   
        for j = 1:102
               UvelPeriodicLiang((i-1)*102 + j) =  Liang(i,j,4);
                             
        end 
 end  
 for i = 1:41   
        for j = 1:102
               VvelPeriodicLiang((i-1)*102 + j) =  Liang(i,j,5);
                             
        end 
 end 
for i = 1:41   
        for j = 1:102
               WvelPeriodicLiang((i-1)*102 + j) =  Liang(i,j,6);
                             
        end 
end  
for i = 1:41   
        for j = 1:102
               LevelSetPeriodicLiang((i-1)*102 + j) =  Liang(i,j,7);
                             
        end 
end 
UvelFinal = UvelPeriodicLiang'; 
VvelFinal = VvelPeriodicLiang'; 
WvelFinal = WvelPeriodicLiang'; 
LevelSetPeriodicLiangFinal = LevelSetPeriodicLiang';

for i = 1:4182
    VelocityFinal(i,1)=UvelFinal(i);
    VelocityFinal(i,2)=VvelFinal(i);
    VelocityFinal(i,3)=WvelFinal(i);
    
    VelocityWithLevelSetFinal(i,1)=UvelFinal(i);
    VelocityWithLevelSetFinal(i,2)=VvelFinal(i);
    VelocityWithLevelSetFinal(i,3)=WvelFinal(i);
    VelocityWithLevelSetFinal(i,4)=LevelSetPeriodicLiangFinal(i);
end