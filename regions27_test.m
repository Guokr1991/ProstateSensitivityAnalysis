
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Method 1: Load data
Pathdata=importdata('/luscinia/ProstateStudy/invivo/Sensitivity_Prostate_Project/Histopathology_27_Regions.csv'); %1a. loads as structure with text 
Pathdata=Pathdata.data; %1b. removes text leaving data only
Pathdata(:,1)=[]; 
ARFIdata=importdata('/luscinia/ProstateStudy/invivo/Sensitivity_Prostate_Project/ARFI_27_Regions.csv');
ARFIdata=ARFIdata.data;
ARFIdata(:,1)=[];
%Method 2: Configure data into standardized form [see Import.m for further
%details]
[Pathdatastandardizedarray,ARFIdatastandardizededarray,Volume,Gleason,Elevation]=Import(Pathdata,ARFIdata);
%Method 3: Configure nearest neighbor matrix, essentially if pathology
%lesion in region/column X, then ARFI ROI allowed to be in regions/columns
%set equivalent to pathology lesion region/column X.
B=zeros(7,27);
B(:,1)=[1 2 14 13 7 19 3];  
B(:,2)=[2 1 13 14 4 0 0]; 
B(:,3)=[3 4 15 9 5 1 0];
B(:,4)=[4 16 3 15 2 6 0];
B(:,5)=[5 6 11 17 27 3 18];
B(:,6)=[6 5 18 17 4 0 0];
B(:,7)=[7 1 8 19 20 13 9];
B(:,8)=[8 7 20 10 19 0 0]
B(:,9)=[9 10 3 21 19 11 0];
B(:,10)=[10 9 22 21 8 12 0];
B(:,11)=[11 12 5 23 27 24 9];
B(:,12)=[12 11 24 23 10 0 0];
B(:,13)=[13 25 14 19 15 0 0];
B(:,14)=[14 13 25 2 16 0 0];
B(:,15)=[15 17 4 3 26 17 13];
B(:,16)=[16 15 26 4 14 18 0];
B(:,17)=[17 28 18 15 5 0 0];
B(:,18)=[18 17 27 6 16 0 0];
B(:,19)=[19 13 25 20 7 21 0];
B(:,20)=[20 8 19 25 22 0 0];
B(:,21)=[21 22 26 19 23 0 0];
B(:,22)=[22 10 20 24 21 0 0];
B(:,23)=[23 24 27 12 21 0 0];
B(:,24)=[24 23 27 12 22 0 0];
B(:,25)=[25 13 19 20 14 26 0];
B(:,26)=[26 15 21 16 22 25 27];
B(:,27)=[27 18 24 17 23 26 0];
m=0;
%Method 4: Takes region location of pathology lesion for patient X,
%compares all ARFI ROIs against this pathology lesion creating a NX1 matrix
%where N is the number of ARFI ROIs. If the NX1 matrix contains a 1, then
%this number is mapped to a 49X1 matrix, and assigned to the given lesion.
%If the NX1 matrix contains only zeros, a zero is mapped to the 49X1
%matrix, and assigned to the given lesion.
for j=1:29
    
Kappa=0; % 4a.Histopathology column number
Phi=0; % 4b.ARFI column number
Pathnumbers=length(find(Pathdatastandardizedarray(:,1,1)==j)) %4c.number of histopathology rows  that satisfy patient number
ARFInumbers=length(find(ARFIdatastandardizededarray(:,1,1)==j)) %4d.number of ARFI rows that satisfy patient number
Pathrows=Pathdatastandardizedarray((find(Pathdatastandardizedarray(:,1,1)==j)),:,1); %4e.storing path rows that correspond to a given patient
Pathrows(:,1)=[]; %4f.deletes patient identifier
ARFIrows=ARFIdatastandardizededarray((find(ARFIdatastandardizededarray(:,1,1)==j)),:,1);
ARFIrows(:,1)=[];
%4f.for loop identifies column numbers (regions) on each path row where cancer
%is located
for i=1:Pathnumbers
Kappa(i)=find(Pathrows(i,:)>0)
end
%4g.for loop identifies column numbers (regions) on each ARFI row where cancer
%is located
for i=1:ARFInumbers
 if isempty(find(ARFIrows(i,:)>0)) %case 1: no call
Phi(i)=-1
   else
Phi(i)=find(ARFIrows(i,:)>0) %case 2: >0 calls
   end
end
%4h. heart of the code: for each histopathology lesion k in Patient j,
%length(Phi) ARFI lesions are compared against it
Mu(j)=length(Kappa); %4i.for each patient number j, Mu identifies number of histopathology lesions that must be
%iterated through
for k=1:Mu(j) 
    m=0; %4j. resets counter when building lesion k specific array
    binary=0; %4k. resets lesion k specific array, note binary=lesion k specific array
   for n=1:length(Phi)
     m=m+1
   binary(m)=(any(B(:,Kappa(k))==Phi(n))); %4l.lesion k specific array
   end
 [Volumerow]=find(binary(:)==1); %4m. searches if 1 contained in binary lesion k specific array
 %4n: based on finding in Volumerow, maps to truth matrix, basically a
 %truth table
 if isempty(Volumerow)==0 & j==1
truth(k,:)=[1 j];
    elseif isempty(Volumerow)==0 & j>1
c=sum(Mu(1:j-1));
truth(c+k,:)=[1 j];
 elseif isempty(Volumerow) & j==1
 truth(k,:)=[0 j];
elseif isempty(Volumerow) & j>1
c=sum(Mu(1:j-1));
truth(c+k,:)=[0 j];
 end
end
end

%Method 5: Filtering data for lesion sub-sets
Overallinformation=zeros(49,5);
Overallinformation(:,:,:,:,:)=[truth(:,1) truth(:,2) Volume' Gleason' Elevation'];
[Posteriorrow Posteriorcolumn]=find(Overallinformation(:,5)==2);
Posteriorrowcolumn=sortrows([Posteriorrow Posteriorcolumn])
Posteriorinformation(:,:,:,:,:)=Overallinformation(Posteriorrow,1:5);
[Anteriorrow Anteriorcolumn]=find(Overallinformation(:,5)==1);
Anteriorrowcolumn=sortrows([Anteriorrow Anteriorcolumn])
Anteriorinformation(:,:,:,:,:)=Overallinformation(Anteriorrow,1:5);
[CSrow CScolumn]=find(Overallinformation(:,3)>500 | Overallinformation(:,4)>6)
CSinformation(:,:,:,:,:)=Overallinformation(CSrow,1:5);
[CINSrow CINScolumn]=find(Overallinformation(:,3)<500 & Overallinformation(:,4)==6 | Overallinformation(:,3)<500 & Overallinformation(:,4)==5);
CINSinformation(:,:,:,:,:)=Overallinformation(CINSrow,1:5);
[CSPosteriorrow CSPosteriorcolumn]=find(Posteriorinformation(:,3)>500 | Posteriorinformation(:,4)>6)
CSPosteriorinformation(:,:,:,:,:)=Posteriorinformation(CSPosteriorrow,1:5);
[CSAnteriorrow CSAnteriorcolumn]=find(Anteriorinformation(:,3)>500 | Anteriorinformation(:,4)>6)
CSAnteriorinformation(:,:,:,:,:)=Anteriorinformation(CSAnteriorrow,1:5);
[CINSPosteriorrow CINSPosteriorcolumn]=find(Posteriorinformation(:,3)<500 & Posteriorinformation(:,4)==6 | Posteriorinformation(:,3)<500 & Posteriorinformation(:,4)==5)
CINSPosteriorinformation(:,:,:,:,:)=Posteriorinformation(CINSPosteriorrow,1:5);
[CINSAnteriorrow CINSAnteriorcolumn]=find(Anteriorinformation(:,3)<500 & Anteriorinformation(:,4)==6 | Anteriorinformation(:,3)<500 & Anteriorinformation(:,4)==5)
CINSAnteriorinformation(:,:,:,:,:)=Overallinformation(CINSAnteriorrow,1:5);
[B I J]=unique(Overallinformation(:,2),'first');
Indexinformation(:,:,:,:,:)=Overallinformation(I,1:5);
[row column]=find(Indexinformation(:,5)==2)
Indexposteriorinformation=Indexinformation(row,1:5);
[row column]=find(Indexinformation(:,5)==1)
Indexanteriorinformation=Indexinformation(row, 1:5);

%Method 6: Sensitivity calculations for lesion subsets, note: Number
%_Identified is the lesion number correctly identified by ARFI, Total_Number is the
%total number of lesions in each subset, Sensitivity was calculated by
%dividing Number_Identified by Total_Number
Number_Identified(:,:,:,:,:)=[sum(Overallinformation(:,1)) sum(Posteriorinformation(:,1)) sum(Anteriorinformation(:,1)); sum(CSinformation(:,1)) sum(CSPosteriorinformation(:,1)) sum(CSAnteriorinformation(:,1));  sum(CINSinformation(:,1))...
    sum(CINSPosteriorinformation(:,1)) sum(CINSAnteriorinformation(:,1)); sum(Indexinformation(:,1))  sum(Indexposteriorinformation(:,1)) sum(Indexanteriorinformation(:,1))]
Total_Number(:,:,:,:,:)=[size(Overallinformation,1) size(Posteriorinformation,1) size(Anteriorinformation,1); size(CSinformation,1)...
    size(CSPosteriorinformation,1) size(CSAnteriorinformation,1); size(CINSinformation,1) size(CINSPosteriorinformation,1)...
   size(CINSAnteriorinformation,1); size(Indexinformation,1) size(Indexposteriorinformation,1) size(Indexanteriorinformation,1)]
Sensitivity(:,:,:,:,:)=Number_Identified./Total_Number;


%Method 7: Graphical and Text Data Display
figure(1)
hBars=bar(Sensitivity);
 Labels = {'Overall', 'Clinically Significant', 'Clinically Insignificant', 'Index Lesions'};
   set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
   legend('Overall Lesions', 'Posterior Lesions', 'Anterior Lesions')
ylabel('Sensitivity')
title('ARFI Sensitivity')


A=char('Overal1', 'Overall Posterior', 'Overall Anterior',  'CS Overall', 'CS Posterior', 'CS Anterior', 'CINS Overall', 'CINS Posterior', 'CINS Anterior', 'Index Overall', 'Index Posterior', 'Index Anterior');
Number_Identified=[Number_Identified(1,:) Number_Identified(2,:) Number_Identified(3,:) Number_Identified(4,:)];
Total_Number=[Total_Number(1,:) Total_Number(2,:) Total_Number(3,:) Total_Number(4,:)];
Sensitivity=[Sensitivity(1,:) Sensitivity(2,:) Sensitivity(3,:) Sensitivity(4,:)]'

fprintf('Sensitivity Values\n\n')
fprintf('Lesion Types\t   Sensitivity\tTP\tTotal Number\n')
for i=1:12
fprintf('%s\t%s\t%d\t%d\n',A(i,:), Sensitivity(i), Number_Identified(i), Total_Number(i))
end
