function [sesInfo] = icatb_dataReductionSetup(sesInfo, optional)
% icatb_dataReductionSetup(sesInfo)
%--------------------------------------------------------------------------
% CREATED BY: Eric Egolf
% LAST MODIFIED: 12-29-03
% ABOUT: Takes data reduction parameters input by user and turns into an
% easier to understand format. It also checks to make sure parameters are
% valid. Used by parameter initialization
%
%
% #########################################################################
%
% USAGE:
%
% -icatb_dataReductionSetup(sesInfo)
%   INFO: sets up variables that make more sense in the context of data
%   reduction than, than the user input variables
%   PARAMETERS:
%       sesInfo = structure containing all parameters nessisary for group
%       ica analysis.
%
% #############################################################
%
% LICENSING:
%
%------------------------------------------------------

modalityType = icatb_get_modality;

if strcmpi(modalityType, 'eeg')
    time_points_type = 'EEG electrodes';
elseif strcmpi(modalityType, 'smri')
    time_points_type = 'Subjects';
else
    time_points_type = 'BOLD timepoints';
end

if ~exist('optional', 'var')
    optional = 'print';
end

% Print information
if strcmpi(optional, 'print')
    disp(' ');
    disp('------------------------------------------------------------------------------------');
    disp('GETTING DATA REDUCTION PARAMETERS----------------------');
    disp('------------------------------------------------------------------------------------');
end

numOfPC1 = sesInfo.userInput.numOfPC1;
numOfPC2 = sesInfo.userInput.numOfPC2;
numOfPC3 = 0;
%numOfPC3 = sesInfo.userInput.numOfPC3;
numOfGroupsStart = sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess; %Change
numOfGroups1 = sesInfo.userInput.numOfGroups1;
numOfGroups2 = sesInfo.userInput.numOfGroups2;
numOfGroups3 = sesInfo.userInput.numOfGroups3;

if(sesInfo.numReductionSteps >0)
    %--SET UP VALUES FOR FIRST REDUCTION STEP
    %---------------------------------------------
    
    %number of groups before and after concatenation
    sesInfo.reduction(1).numOfGroupsBeforeCAT=numOfGroupsStart;
    sesInfo.reduction(1).numOfGroupsAfterCAT=numOfGroups1;
    %number of principal componentes before and after principal component
    %analysis
    % taking care of the different time points
    sesInfo.reduction(1).numOfPCBeforeCAT = sesInfo.diffTimePoints; %sesInfo.numOfScans;
    sesInfo.reduction(1).numOfPCAfterReduction=numOfPC1;
    %number of groups in each newly concatenated group
    sesInfo.reduction(1).numOfPrevGroupsInEachNewGroupAfterCAT([1:sesInfo.reduction(1).numOfGroupsAfterCAT])=1;
    %number of principal components in each group
    sesInfo.reduction(1).numOfPCInEachGroupAfterCAT([1:sesInfo.reduction(1).numOfGroupsAfterCAT])=sesInfo.reduction(1).numOfPrevGroupsInEachNewGroupAfterCAT.*sesInfo.reduction(1).numOfPCBeforeCAT;
    
    
    
    
    
    %--SET UP VALUES FOR SECOND REDUCTION STEP
    %-------------------------------------------------------
    %check if oddman out case
    groupsBefore = sesInfo.reduction(1).numOfGroupsAfterCAT;
    newGroups = numOfGroups2;
    isOddManOut = mod(groupsBefore,newGroups);
    if(~isOddManOut)
        
        sesInfo.reduction(2).numOfGroupsBeforeCAT=sesInfo.reduction(1).numOfGroupsAfterCAT;
        sesInfo.reduction(2).numOfGroupsAfterCAT=numOfGroups2;
        
        sesInfo.reduction(2).numOfPCBeforeCAT=numOfPC1;
        sesInfo.reduction(2).numOfPCAfterReduction=numOfPC2;
        
        sesInfo.reduction(2).numOfPrevGroupsInEachNewGroupAfterCAT([1:sesInfo.reduction(2).numOfGroupsAfterCAT])=...
            sesInfo.reduction(2).numOfGroupsBeforeCAT/sesInfo.reduction(2).numOfGroupsAfterCAT;
        
        sesInfo.reduction(2).numOfPCInEachGroupAfterCAT([1:sesInfo.reduction(2).numOfGroupsAfterCAT])=...
            sesInfo.reduction(2).numOfPCBeforeCAT*sesInfo.reduction(2).numOfPrevGroupsInEachNewGroupAfterCAT(1);
        
    else
        
        sesInfo.reduction(2).numOfGroupsBeforeCAT=numOfGroups1;
        sesInfo.reduction(2).numOfGroupsAfterCAT=numOfGroups2-1;
        
        sesInfo.reduction(2).numOfPCBeforeCAT=numOfPC1;
        sesInfo.reduction(2).numOfPCAfterReduction=numOfPC2;
        %----------------------------------------------------------------------
        %Oddman out case
        %----------------------------------------------------------------------
        totalPrevGroupsPutInNewGroups=0;
        for i=1:newGroups
            % Last group from pca does not fit evenly into new
            % group
            if(i~=newGroups)
                
                sesInfo.reduction(2).numOfPrevGroupsInEachNewGroupAfterCAT(i)=...
                    floor(sesInfo.reduction(2).numOfGroupsBeforeCAT/sesInfo.reduction(2).numOfGroupsAfterCAT);
                
                sesInfo.reduction(2).numOfPCInEachGroupAfterCAT(i)=...
                    sesInfo.reduction(1).numOfPCAfterReduction*sesInfo.reduction(2).numOfPrevGroupsInEachNewGroupAfterCAT(i);
                
            else
                
                % Put oddman out in previous group
                numLeftOut = sesInfo.reduction(2).numOfGroupsBeforeCAT - totalPrevGroupsPutInNewGroups;
                
                sesInfo.reduction(2).numOfPrevGroupsInEachNewGroupAfterCAT(i-1)=sesInfo.reduction(2).numOfPrevGroupsInEachNewGroupAfterCAT(i-1)+numLeftOut;
                
                sesInfo.reduction(2).numOfPCInEachGroupAfterCAT(i-1)=sesInfo.reduction(2).numOfPCInEachGroupAfterCAT(i-1)+numLeftOut*sesInfo.reduction(2).numOfPCBeforeCAT;
                
                % Display feedback for user that their
                % specifications are being changed
                string=['    NOTE: After reduction step ',num2str(2) ,' can not make ',num2str(newGroups),' full groups from ',num2str(sesInfo.reduction(2).numOfGroupsBeforeCAT)];
                string2=['    Alternatively making ', num2str(newGroups-1),' groups from ',num2str(sesInfo.reduction(2).numOfGroupsBeforeCAT)];
                if strcmpi(optional, 'print')
                    disp(string);
                    disp(string2);
                end
            end
            
            
            %In case last run of oddman out, don't calculate total
            if(i~=newGroups)
                totalPrevGroupsPutInNewGroups=totalPrevGroupsPutInNewGroups+sesInfo.reduction(2).numOfPrevGroupsInEachNewGroupAfterCAT(i);
            end
            
        end
        
    end
    
    
    
    
    
    
    %--SET UP VALUES FOR THIRD REDUCTION STEP
    %-------------------------------------------------------
    %check if oddman out case
    groupsBefore = sesInfo.reduction(2).numOfGroupsAfterCAT;
    newGroups = numOfGroups3;
    isOddManOut = mod(groupsBefore,newGroups);
    if(~isOddManOut)
        
        sesInfo.reduction(3).numOfGroupsBeforeCAT=sesInfo.reduction(2).numOfGroupsAfterCAT;
        sesInfo.reduction(3).numOfGroupsAfterCAT=numOfGroups3;
        
        sesInfo.reduction(3).numOfPCBeforeCAT=numOfPC2;
        sesInfo.reduction(3).numOfPCAfterReduction=numOfPC3;
        
        sesInfo.reduction(3).numOfPrevGroupsInEachNewGroupAfterCAT([1:sesInfo.reduction(3).numOfGroupsAfterCAT])=...
            sesInfo.reduction(3).numOfGroupsBeforeCAT/sesInfo.reduction(3).numOfGroupsAfterCAT;
        
        sesInfo.reduction(3).numOfPCInEachGroupAfterCAT([1:sesInfo.reduction(3).numOfGroupsAfterCAT])=...
            sesInfo.reduction(3).numOfPCBeforeCAT*sesInfo.reduction(3).numOfPrevGroupsInEachNewGroupAfterCAT(1);
        
    else
        
        sesInfo.reduction(3).numOfGroupsBeforeCAT=sesInfo.reduction(2).numOfGroupsAfterCAT;
        sesInfo.reduction(3).numOfGroupsAfterCAT=numOfGroups3-1;
        
        sesInfo.reduction(3).numOfPCBeforeCAT=numOfPC2;
        sesInfo.reduction(3).numOfPCAfterReduction=numOfPC3;
        %----------------------------------------------------------------------
        %Oddman out case
        %----------------------------------------------------------------------
        totalPrevGroupsPutInNewGroups=0;
        for i=1:newGroups
            % Last group from pca does not fit evenly into new
            % group
            if(i~=newGroups)
                
                sesInfo.reduction(3).numOfPrevGroupsInEachNewGroupAfterCAT(i)=...
                    floor(sesInfo.reduction(3).numOfGroupsBeforeCAT/sesInfo.reduction(3).numOfGroupsAfterCAT);
                
                sesInfo.reduction(3).numOfPCInEachGroupAfterCAT(i)=...
                    sesInfo.reduction(2).numOfPCAfterReduction*sesInfo.reduction(3).numOfPrevGroupsInEachNewGroupAfterCAT(i);
                
            else
                
                % Put oddman out in previous group
                numLeftOut = sesInfo.reduction(3).numOfGroupsBeforeCAT - totalPrevGroupsPutInNewGroups;
                
                sesInfo.reduction(3).numOfPrevGroupsInEachNewGroupAfterCAT(i-1)=sesInfo.reduction(3).numOfPrevGroupsInEachNewGroupAfterCAT(i-1)+numLeftOut;
                
                sesInfo.reduction(3).numOfPCInEachGroupAfterCAT(i-1)=sesInfo.reduction(3).numOfPCInEachGroupAfterCAT(i-1)+numLeftOut*sesInfo.reduction(3).numOfPCBeforeCAT;
                
                % Display feedback for user that their
                % specifications are being changed
                string=['    NOTE: After reduction step ',num2str(3) ,' can not make ',num2str(newGroups),' full groups from ',num2str(sesInfo.reduction(3).numOfGroupsBeforeCAT)];
                string2=['    Alternatively making ', num2str(newGroups-1),' groups from ',num2str(sesInfo.reduction(3).numOfGroupsBeforeCAT)];
                if strcmpi(optional, 'print')
                    disp(string);
                    disp(string2);
                end
            end
            
            
            %In case last run of oddman out, don't calculate total
            if(i~=newGroups)
                totalPrevGroupsPutInNewGroups=totalPrevGroupsPutInNewGroups+sesInfo.reduction(2).numOfPrevGroupsInEachNewGroupAfterCAT(i);
            end
            
        end
        
    end
    
    
    
    
    
    for j=1:sesInfo.numReductionSteps
        string =[' Reduction step ',num2str(j),' starts with ',num2str(sesInfo.reduction(j).numOfGroupsBeforeCAT),' groups and gets reduced to ',num2str(sesInfo.reduction(j).numOfGroupsAfterCAT), ' groups'];
        if strcmpi(optional, 'print')
            disp(string);
        end
        
        for i =1:sesInfo.reduction(j).numOfGroupsAfterCAT
            string=['    -New Group #', num2str(i), ': ',num2str(sesInfo.reduction(j).numOfPrevGroupsInEachNewGroupAfterCAT(i)),' groups will be concatenated to form the new group.'];
            if j==1
                string2=['     Each of the to be concatenated groups is made up of ', num2str(sesInfo.reduction(j).numOfPCBeforeCAT(i)), ' ', time_points_type];
                string3=['     The new group will have a total of ', num2str(sesInfo.reduction(j).numOfPCBeforeCAT(i)*sesInfo.reduction(j).numOfPrevGroupsInEachNewGroupAfterCAT(i)),' stacked ', time_points_type];
            else
                string2=['     Each of the to be concatenated groups is made up of ', num2str(sesInfo.reduction(j).numOfPCBeforeCAT),' principal components'];
                string3=['     The new group will have a total of ', num2str(sesInfo.reduction(j).numOfPCBeforeCAT*sesInfo.reduction(j).numOfPrevGroupsInEachNewGroupAfterCAT(i)),' stacked principal components'];
            end
            
            string4=['     This group will then be reduced to ',num2str(sesInfo.reduction(j).numOfPCAfterReduction), ' principal components'];
            if strcmpi(optional, 'print')
                disp(string);
                disp(string2);
                disp(string3);
                disp(string4);
            end
        end
    end
else
    if strcmpi(optional, 'print')
        disp('--No Data Reduction Steps');
    end
end

if strcmpi(optional, 'print')
    disp('------------------------------------------------------------------------------------');
    disp('END GETTING DATA REDUCTION PARAMETERS----------------------');
    disp('------------------------------------------------------------------------------------');
    disp(' ');
end