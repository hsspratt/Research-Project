function [fitresult, gof] = createFit(LambdaStraightLine, StraightLine)
%CREATEFIT(LAMBDASTRAIGHTLINE,STRAIGHTLINE)
%  Create a fit.
%
%  Data for 'Straight Line Best Fit' fit:
%      X Input : LambdaStraightLine
%      Y Output: StraightLine
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 04-Jan-2022 23:31:20

%% Fit: 'Straight Line Best Fit'.
[xData, yData] = prepareCurveData( LambdaStraightLine, StraightLine );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Create a figure for the plots.
figure( 'Name', 'Straight Line Best Fit' );

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult, xData, yData, 'predobs', 0.9 );
legend( h, 'StraightLine vs. LambdaStraightLine', 'Straight Line Best Fit', 'Lower bounds (Straight Line Best Fit)', 'Upper bounds (Straight Line Best Fit)', 'Location', 'NorthWest', 'Interpreter', 'none' );
% Label axes
xlabel( 'LambdaStraightLine', 'Interpreter', 'none' );
ylabel( 'StraightLine', 'Interpreter', 'none' );
grid on

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult, xData, yData, 'residuals' );
legend( h, 'Straight Line Best Fit - residuals', 'Zero Line', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'LambdaStraightLine', 'Interpreter', 'none' );
ylabel( 'StraightLine', 'Interpreter', 'none' );
grid on


