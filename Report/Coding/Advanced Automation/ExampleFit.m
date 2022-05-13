[fitresult, gof] = createFit1(energy_ev, alpha_1)

plot( fitresult, energy_ev, alpha_1, 'predobs', 0.99 )

function [fitresult, gof] = createFit1(energy_ev, alpha_1)
%CREATEFIT1(ENERGY_EV,ALPHA_1)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : energy_ev
%      Y Output: alpha_1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 13-May-2022 15:03:57


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( energy_ev, alpha_1 );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
%figure( 'Name', 'untitled fit 1' );
%h = plot( fitresult, xData, yData, 'predobs', 0.99 );
%legend( h, 'alpha_1 vs. energy_ev', 'untitled fit 1', 'Lower bounds (untitled fit 1)', 'Upper bounds (untitled fit 1)', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
%xlabel( 'energy_ev', 'Interpreter', 'none' );
%ylabel( 'alpha_1', 'Interpreter', 'none' );
%grid on
end
