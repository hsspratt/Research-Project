%% Optical Anlysis Functions

classdef OpticalAnalysisFunctions
    methods(Static)
        % Refraction True Value

        function RefractionTrueValue = nearestRefraction(LambdaStore, Lambda, RefractionIndex)
        
        lambda_d = Lambda-LambdaStore;
        
        if isempty(find(lambda_d==0)) == 0
            lowest=find(lambda_d==0);
        elseif isempty(find(lambda_d==0)) == 1
            lowest=max(find(lambda_d>0));
        end
        
        highest=lowest+1;
        
        n_d = RefractionIndex(highest)- RefractionIndex(lowest);
        
        if n_d ~= 0
            Lambda_Percentage = lambda_d(lowest)./(LambdaStore(highest)-LambdaStore(lowest));
        elseif n_d == 0
            Lambda_Percentage = 0;
        end
        
        if n_d > 0
            RefractionTrueValue = RefractionIndex(lowest) - (n_d).* Lambda_Percentage;
        elseif n_d < 0
            RefractionTrueValue = RefractionIndex(lowest) + (n_d).* Lambda_Percentage;
        elseif n_d == 0
            RefractionTrueValue = RefractionIndex(lowest);
        else 
            disp('Error Alert')
        end
        end
        
        % Straight Line Values - Could be multiple sections

        function StraightLine = DetectStraightLine(xdata, ydata)
        
        gradient(1) = ydata(1)./(xdata(1));
        for i=2:size(xdata,2)-1
            gradient(i) = (ydata(i+1) - ydata(i))./(xdata(i+1) - xdata(i));
        end

        gradient = smoothdata(gradient);
        StraightLine = [];
        for i = 1:size(gradient,2)-1
            section = gradient(i:i+1);
            p12 = section(1);
            p23 = section(2);
            per = p23*0.05;
        
            if (p23-per <= p12) && (p12 <= p23+per)
               StraightLine(i) = section(1);
            end
        end

        if isempty(StraightLine) == 1
            disp('There is no clear straightline')
        end
        end

        % Longest Consecutive straight line and index - UPDATED

        function [LongestConsecutive, longestLengthIndex] = LongestConsecutive(X)
        
        longestLength = 0;
        longestLengthIndex = [];
        
        currLength = 1;
        
        for i = 1:size(X,2)-1
            if X(i+1) == X(i) + 1
                currLength = currLength + 1;
                longestLength = max(longestLength, currLength);
            else
                currLength = 1;
                
            end
        LongestConsecutive = max(longestLength, currLength);
        if currLength == longestLength
            longestLengthIndex = i+1;
        end
        
        end
        end
        
        function [xdata, ydata] = DetectLongestStarightLine(xdata, ydata)

            StraightLine = OpticalAnalysisFunctions.DetectStraightLine(xdata, ydata); % detects staright line sections

            nonzero = find(StraightLine ~= 0); % finds nonstraight points from array
            
            %x_data_excess = x_data(nonzero); % removes nonstraight points from array
            %y_data_excess = y_data(nonzero); % removes nonstraight points from array

            [length, index] = OpticalAnalysisFunctions.LongestConsecutive(nonzero); % finds longest section + index
            nonzero_corrected = nonzero(index-length+1:index); % corrects it by a factor of + 1
            
            xdata = xdata(nonzero_corrected);
            ydata = ydata(nonzero_corrected);


        end
    end
end





% %% Optical Anlysis Functions
% 
% % Refraction True Value
% 
% function RefractionTrueValue = nearestRefraction(LambdaStore, Lambda, RefractionIndex)
% 
% lambda_d = Lambda-LambdaStore;
% 
% if isempty(find(lambda_d==0)) == 0
%     lowest=find(lambda_d==0);
% elseif isempty(find(lambda_d==0)) == 1
%     lowest=max(find(lambda_d>0));
% end
% 
% highest=lowest+1;
% 
% n_d = RefractionIndex(highest)- RefractionIndex(lowest);
% 
% if n_d ~= 0
%     Lambda_Percentage = lambda_d(lowest)./(LambdaStore(highest)-LambdaStore(lowest));
% elseif n_d == 0
%     Lambda_Percentage = 0;
% end
% 
% if n_d > 0
%     RefractionTrueValue = RefractionIndex(lowest) - (n_d).* Lambda_Percentage;
% elseif n_d < 0
%     RefractionTrueValue = RefractionIndex(lowest) + (n_d).* Lambda_Percentage;
% elseif n_d == 0
%     RefractionTrueValue = RefractionIndex(lowest);
% else 
%     disp('Error Alert')
% end
% end
% 
% % Straight Line Values - Could be multiple sections
% 
% function StraightLine = DetectStraightLine(data)
% data = smoothdata(data);
% StraightLine = [];
% for i = 1:size(data,2)
%     section = data(i:i+1);
%     p12 = section(i);
%     p23 = section(i+1);
%     per = p23*0.05;
% 
%     if (p23-per <= p12) && (p12 <= p23+per)
%        StraightLine(i) = section(i);
%     end
% end
% end
% 
% % Longest Consecutive straight line and index - UPDATED
% 
% function [LongestConsecutive, longestLengthIndex] = LongestConsecutive2(X)
% 
% longestLength = 0;
% longestLengthIndex = [];
% 
% currLength = 1;
% 
% for i = 1:size(X,2)-1
%     if X(i+1) == X(i) + 1
%         currLength = currLength + 1;
%         longestLength = max(longestLength, currLength);
%     else
%         currLength = 1;
%         
%     end
% LongestConsecutive = max(longestLength, currLength);
% if currLength == longestLength
%     longestLengthIndex = i+1;
% end
% 
% end
% end

