%% Optical Anlysis Functions

classdef OpticalAnalysisFunctions
    methods(Static)
        % y True Value

        function y_true_value = nearestValue(x_Store, x, y_Store)
        
        x_diffs = x-x_Store;
        
        if isempty(find(x_diffs==0)) == 0
            % zero_diff=find(x_diffs==0);
            % lowest = (min(zero_diff)+max(zero_diff))/length(zero_diff)
            lowest = min(find(x_diffs==0))
        elseif isempty(find(x_diffs==0)) == 1
            lowest=max(find(x_diffs>0));
        end
        
        highest=lowest+1;
        
        y_diff = y_Store(highest)- y_Store(lowest);
        
        if y_diff ~= 0
            x_Percentage = x_diffs(lowest)./(x_Store(highest)-x_Store(lowest));
        elseif y_diff == 0
            x_Percentage = 0;
        end
        
        if y_diff ~= 0
            y_true_value = y_Store(lowest) + (y_diff).* x_Percentage;
        elseif y_diff == 0
            y_true_value = y_Store(lowest);
        else 
            disp('Error Alert')
        end
        end
        
        % Straight Line Values - Could be multiple sections

        function StraightLine = DetectStraightLine(xdata, ydata, err)
        
        xdata = smooth(xdata);
        ydata = smooth(ydata);

        gradient = zeros(max(size(xdata))-1, 1);
        gradient(1) = ydata(1)./(xdata(1));
        for i=1:max(size(xdata))-1
            gradient(i) = (ydata(i+1) - ydata(i))./(xdata(i+1) - xdata(i));
        end

        gradient = smoothdata(gradient);
        StraightLine = [];
        for i = 1:size(gradient,1)-1
            section = gradient(i:i+1);
            if abs(gradient) > (max(ydata)-min(ydata))/((max(xdata)-min(xdata)))
                p12 = section(1);
                p23 = section(2);
                % per = p23*err;
                rolling_average = mean(gradient(1:i+1));
                per = rolling_average*err;
    
                if (p23-per <= p12) && (p12 <= p23+per)
                   StraightLine(i) = section(1);
                end
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
        
        function [xdata, ydata] = DetectLongestStarightLine(xdata, ydata, err)

            StraightLine = OpticalAnalysisFunctions.DetectStraightLine(xdata, ydata, err); % detects staright line sections

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
            per = max(ydata)*(10/100);
        
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

        function StraightLine = DetectStraightLine1(xdata, ydata, err)
        
        gradient(1) = ydata(1)./(xdata(1));
        for i=1:max(size(xdata))-1
            gradient(i) = (ydata(i+1) - ydata(i))./(xdata(i+1) - xdata(i));
        end

        gradient = smoothdata(gradient);
        min_gradient = (ydata(end) - ydata(1))/(xdata(end) - xdata(1));
        
        for i=1:max(size(gradient))
            if gradient(i)<min_gradient
                gradient(i)=[]
            end
        end
 
        StraightLine = [];
        for i = 1:size(gradient,2)-1
            section = gradient(i:i+1);
            p12 = section(1);
            p23 = section(2);
            %per = p23*err;
            rolling_average = mean(gradient(1:i+1))
            per = rolling_average*err;

            if (p23-per <= p12) && (p12 <= p23+per)
               StraightLine(i) = section(1);
            end
        end

        if isempty(StraightLine) == 1
            disp('There is no clear straightline')
        end
        end

        function RefractiveIndex = CalculateRefractiveIndex(wavelength, type)

            if strcmp(type,'GaAs') == 1
                RefractiveIndex = csvimport('Papatryfonoset-2021-0.260-1.88-Ga_As.csv');
            elseif strcmp(type,'GaP')
                disp('Need to connect csv file!!')
            end

            RefractiveIndex(1,:)  = [];  % removing column titles
            
            RefractiveIndexInfo  = zeros(max(size(RefractiveIndex)),min(size(RefractiveIndex)));
            
            for c=1:min(size(RefractiveIndex))
                for r=1:max(size(RefractiveIndex))
                    RefractiveIndexInfo(r,c)  = RefractiveIndex{r,c};
                end
            end
            
            L  = RefractiveIndexInfo(:,1).*1000;
            N  = RefractiveIndexInfo(:,2);
            K  = RefractiveIndexInfo(:,3);
            
            n = zeros(size(wavelength));
            k = zeros(size(wavelength));
            R = zeros(size(wavelength));
            
            for i=1:max(size(wavelength))
                % Approximated value of R calculated - Not calculating R manually with this method
                n(i) = OpticalAnalysisFunctions.nearestValue(L, wavelength(i), N);
                k(i) = OpticalAnalysisFunctions.nearestValue(L, wavelength(i), K);
                R(i) = ((n(i)-1)+k(i).^2)/((n(i)+1)+k(i).^2);
            end

            RefractiveIndex = R;

            figure( 'Name', 'Refractive Index');

            plot(L,N) % real refractive index - n
            hold on
            plot(L,K) % complex refracive index - ik
            hold on
            plot(wavelength(1:50:end),n(1:50:end),'*') % real refractive index - n
            hold on
            plot(wavelength(1:50:end),k(1:50:end),'*') % complex refracive index - ik
            
            title('Refractive Index $200 - 830 nm$','Interpreter','latex');
            legend('Real refractive index','Complex refractive index', ...
                'Interpolated Real refractive','Interpolated Complex refractive','Interpreter','latex')
            
            % legand needs changing
            
            xlabel( 'Wavelength $(nm)$', 'Interpreter', 'latex' );
            ylabel( 'Refractive Index ', 'Interpreter', 'Latex' );

        end

        function NewWavelengths = WavelengthsSystematicCorrection(wavelengths) 
        m = -0.0091;
        c = 18.1;
        x = wavelengths;
        y = x + (m*x + c);
        % wavelengths = y;
        NewWavelengths = y;
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

