


function [Pathdatastandardizedarray,ARFIdatastandardizedarray,Volume,Gleason,Elevation]=Import(Pathdata,ARFIdata)
%Method: 1. Load data
Pathdatarows=size(Pathdata,1); 
%Method 2: Build standardized array, will be used in later code to put
%original Pathdata array in standard form
Pathdatastandardizedarray=zeros(Pathdatarows,28,4); 
%Method 3: construct 49X28X4 matrix from original array
for i=1:28 
    %3a: sets patient identifiers as first column in each dim of
    %standardized matrix
    if i==1
    for j=1:4
Pathdatastandardizedarray(:,1,j)=Pathdata(:,1); 
    end
 %3b: copies volume+gleason score data to appropriate dimensions
    elseif i>1
Pathdatastandardizedarray(:,i,1)=Pathdata(:,2.*(i-1));  
 Pathdatastandardizedarray(:,i,2)=Pathdata(:,2.*(i-1));
  Pathdatastandardizedarray(:,i,3)=Pathdata(:,2.*i-1); 
 %Method 4: uses volume and gleason score data to construct binary data in
 %dim 1 and posterior/mid-anterior data in dim 4
 %4a. conversion of volumetric data to binary data
   [pathrowoccupied pathcolumnoccupied]=find(Pathdatastandardizedarray(:,i,2)>0);
   if Pathdatastandardizedarray(pathrowoccupied, i,1)>0
Pathdatastandardizedarray(pathrowoccupied,i,1)=1;
   end
   %4b: based on tumor region, assigns different numerical values to
   %posterior colums vs. mid/anterior columns [indexed by the for loop]
    if i<=13
Pathdatastandardizedarray(pathrowoccupied,i,4)=2;
   elseif i>13
Pathdatastandardizedarray(pathrowoccupied,i,4)=1;
    end
    end
end
%Method 5: collect volume/gleason/elevation data in individualized rows for
%easy viewing
  [row column]=find(Pathdatastandardizedarray(:,2:28,2)>0); %5a. finds occupied cells
column=column+ones(49,1); %5b: must take into account 2:28
rowcolumn=[row column]; 
rowcolumn=sortrows(rowcolumn); %5c: sort in row-based numerical order i.e. patient identifier order
%5d: for loop pulls occupied cells out in patient identifier order to
%create individual volume/gleason/elevation columns
%question:is for loop necessary?
for i=1:49
Volume(i)=Pathdatastandardizedarray(rowcolumn(i,1),rowcolumn(i,2),2); 
Gleason(i)=Pathdatastandardizedarray(rowcolumn(i,1),rowcolumn(i,2),3);
Elevation(i)=Pathdatastandardizedarray(rowcolumn(i,1),rowcolumn(i,2),4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Repeat process for ARFI data, modified as only two dimensions
%dim 1: binary data
%dim 2: index of suspicion data
ARFIdatarows=size(ARFIdata,1);
ARFIdatastandardizedarray=zeros(ARFIdatarows, 28,2);

for k=1:28
    if k==1
for j=1:2
ARFIdatastandardizedarray(:,:,j)=ARFIdata(:,:);
end
    elseif k>1
 [ARFIrowoccupied ARFIcolumnoccupied]=find(ARFIdatastandardizedarray(:,k,2)>0);
if ARFIdatastandardizedarray(ARFIrowoccupied,k,1)>0
ARFIdatastandardizedarray(ARFIrowoccupied,k,1)=1;
end
    end
end
