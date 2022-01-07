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
        for i=1:max(size(xdata))-1
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

        function [Transparent_Wavelength, Transparent_T] = TransparentStraightLine(xdata, ydata)
        
        StraightLine = [];
        for i=max(size(ydata)):-1:2
            section = ydata(i:-1:i-1);
            p3 = section(1);
            p2 = section(2);
            per = p3*0.01;
        
            rolling_average = mean(ydata(i:max(size(ydata))));
        
            if (rolling_average-per <= p2) && (p2 <= rolling_average+per)
               StraightLine(i) = ydata(1);
            end
        end
        
        nonzero = find(StraightLine ~= 0);
        [length, index] = OpticalAnalysisFunctions.LongestConsecutive(nonzero);
        nonzero_corrected = nonzero(index-length+1:index);
        Transparent_Wavelength = xdata(nonzero_corrected);
        Transparent_T = ydata(nonzero_corrected);

        end

        function [Excess_Wavelength, Excess_T, wavelength, T, Index] = CutExcessData(xdata, ydata)
        
        StraightLine = [];
        for i=1:1:max(size(ydata))-1
            points = ydata(i:1:i+1);
            p1 = points(1);
            p2 = points(2);
            %per = p1*2;
            per = max(ydata)*(3.5/100);
        
            rolling_average = mean(ydata(1:i+1));
        
            if (rolling_average-per <= p2) && (p2 <= rolling_average+per)
               StraightLine(i) = ydata(1);
            end
        end
        
        nonzero           = find(StraightLine ~= 0);
        [length, index]   = OpticalAnalysisFunctions.LongestConsecutive(nonzero);
        nonzero_corrected = nonzero(index-length+1:index);
        Excess_Wavelength = xdata(nonzero_corrected);
        Excess_T          = ydata(nonzero_corrected);

        Index      = round(max(size(Excess_T))*0.90,0); % Threshold at 80%

        wavelength = xdata(Index:1:end);
        T          = ydata(Index:1:end);

        end

%         function [Excess_Wavelength, Excess_Alpha, wavelength, T, Index] = CutExcessDataAlpha(xdata, ydata)
%         
%         StraightLine = [];
%         for i=max(size(ydata)):-1:2
%             points = ydata(i:-1:i-1);
%             p1 = points(1);
%             p2 = points(2);
%             %per = p1*2;
%             per = sort(ydata);
%             per = mean(per(round(max(size(ydata)),0)))
%         
%             rolling_average = mean(ydata(i:max(size(ydata))));
%         
%             if (rolling_average-per <= p2) && (p2 <= rolling_average+per)
%                StraightLine(i) = ydata(1);
%             end
%         end
%         
%         nonzero           = find(StraightLine ~= 0);
%         [length, index]   = OpticalAnalysisFunctions.LongestConsecutive(nonzero);
%         nonzero_corrected = nonzero(index-length+1:index);
%         Excess_Wavelength = xdata(nonzero_corrected);
%         Excess_Alpha          = ydata(nonzero_corrected);
% 
%         Index      = round(max(size(Excess_Alpha))*0.75,0); % Threshold at 80%
% 
%         wavelength = xdata(Index:1:end);
%         T          = ydata(Index:1:end);
% 
%         end

        function AnalysisData = AutomatedAnalysis()
            
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

