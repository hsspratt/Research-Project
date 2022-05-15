classdef AutomatedAppBackup < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        SigmoidPlot                    matlab.ui.control.CheckBox
        TaucPlot                       matlab.ui.control.CheckBox
        Step                           matlab.ui.control.NumericEditField
        WavelengthStepLabel            matlab.ui.control.Label
        Readings                       matlab.ui.control.NumericEditField
        NumberofMeasurementsLabel      matlab.ui.control.Label
        FastMethod                     matlab.ui.control.CheckBox
        Transmission_eV                matlab.ui.control.CheckBox
        FinalWavelength                matlab.ui.control.NumericEditField
        FinalWavelengthEditFieldLabel  matlab.ui.control.Label
        StartingWavelength             matlab.ui.control.NumericEditField
        StartingWavelengthEditFieldLabel  matlab.ui.control.Label
        RunScriptButton                matlab.ui.control.StateButton
        GraphsPlotted_Title            matlab.ui.control.TextArea
        Parameters_Title               matlab.ui.control.TextArea
        Title_Text                     matlab.ui.control.TextArea
    end

   

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: RunScriptButton
        function RunScriptButtonValueChanged(app, event)
            %value = app.RunScriptButton.Value;
            s_lambda = app.StartingWavelength.Value;
            f_lambda = app.FinalWavelength.Value;
            Readings = app.Readings.Value;
            Step = app.Step.Value;

            Transmission = app.Transmission_eV.Value;
            FastMethod = app.FastMethod.Value;
            Tauc = app.TaucPlot.Value
            Sigmoid = app.SigmoidPlot.Value


            matrix1 = {'Starting Wavelength';'Final Wavelength';
                       'Readings'; 'Steps';
                       'Transmission (eV)'; 'Fast Method';
                        'Tauc Plotting'; 'Sigmoid Fit'};
            matrix2 = [s_lambda;f_lambda;Readings;Step;Transmission;
                       FastMethod;Tauc;Sigmoid];
            fid = fopen( 'Parameters.csv', 'w' );
            for jj = 1 : length( matrix1 )
                fprintf( fid, '%s,%d\n', matrix1{jj}, matrix2(jj) );
            end
            fclose( fid )

            uialert(app.UIFigure,'Upon pressing program will be run', ...
            'Program Information','Icon','info','CloseFcn','uiresume(app.UIFigure)')

            uiwait(app.UIFigure)

            delete(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [1 1 1];
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create Title_Text
            app.Title_Text = uitextarea(app.UIFigure);
            app.Title_Text.FontSize = 40;
            app.Title_Text.Position = [43 375 556 56];
            app.Title_Text.Value = {'Automated Optical Absorption'};

            % Create Parameters_Title
            app.Parameters_Title = uitextarea(app.UIFigure);
            app.Parameters_Title.FontSize = 14;
            app.Parameters_Title.FontWeight = 'bold';
            app.Parameters_Title.Position = [43 332 89 26];
            app.Parameters_Title.Value = {'Parameters'};

            % Create GraphsPlotted_Title
            app.GraphsPlotted_Title = uitextarea(app.UIFigure);
            app.GraphsPlotted_Title.FontSize = 14;
            app.GraphsPlotted_Title.FontWeight = 'bold';
            app.GraphsPlotted_Title.Position = [408 333 112 24];
            app.GraphsPlotted_Title.Value = {'Graphs Plotted'};

            % Create RunScriptButton
            app.RunScriptButton = uibutton(app.UIFigure, 'state');
            app.RunScriptButton.ValueChangedFcn = createCallbackFcn(app, @RunScriptButtonValueChanged, true);
            app.RunScriptButton.Text = 'Run Script';
            app.RunScriptButton.Position = [262 44 100 22];

            % Create StartingWavelengthEditFieldLabel
            app.StartingWavelengthEditFieldLabel = uilabel(app.UIFigure);
            app.StartingWavelengthEditFieldLabel.HorizontalAlignment = 'right';
            app.StartingWavelengthEditFieldLabel.Position = [43 287 114 22];
            app.StartingWavelengthEditFieldLabel.Text = 'Starting Wavelength';

            % Create StartingWavelength
            app.StartingWavelength = uieditfield(app.UIFigure, 'numeric');
            app.StartingWavelength.Position = [202 287 37 22];

            % Create FinalWavelengthEditFieldLabel
            app.FinalWavelengthEditFieldLabel = uilabel(app.UIFigure);
            app.FinalWavelengthEditFieldLabel.HorizontalAlignment = 'right';
            app.FinalWavelengthEditFieldLabel.Position = [43 246 97 22];
            app.FinalWavelengthEditFieldLabel.Text = 'Final Wavelength';

            % Create FinalWavelength
            app.FinalWavelength = uieditfield(app.UIFigure, 'numeric');
            app.FinalWavelength.Position = [202 246 37 22];
            app.FinalWavelength.Value = 1100;

            % Create Transmission_eV
            app.Transmission_eV = uicheckbox(app.UIFigure);
            app.Transmission_eV.Text = '       Transmission (eV)';
            app.Transmission_eV.Position = [408 287 138 22];

            % Create FastMethod
            app.FastMethod = uicheckbox(app.UIFigure);
            app.FastMethod.Text = '       Trivial Band Gap Calculation';
            app.FastMethod.Position = [408 246 198 22];

            % Create NumberofMeasurementsLabel
            app.NumberofMeasurementsLabel = uilabel(app.UIFigure);
            app.NumberofMeasurementsLabel.HorizontalAlignment = 'right';
            app.NumberofMeasurementsLabel.Position = [43 205 145 22];
            app.NumberofMeasurementsLabel.Text = 'Readings per wavelength ';

            % Create Readings
            app.Readings = uieditfield(app.UIFigure, 'numeric');
            app.Readings.Position = [202 205 37 22];
            app.Readings.Value = 4;

            % Create WavelengthStepLabel
            app.WavelengthStepLabel = uilabel(app.UIFigure);
            app.WavelengthStepLabel.HorizontalAlignment = 'right';
            app.WavelengthStepLabel.Position = [43 167 97 22];
            app.WavelengthStepLabel.Text = 'Wavelength Step';

            % Create Step
            app.Step = uieditfield(app.UIFigure, 'numeric');
            app.Step.Position = [202 167 37 22];
            app.Step.Value = 1;

            % Create TaucPlot
            app.TaucPlot = uicheckbox(app.UIFigure);
            app.TaucPlot.Text = '       Tauc Plotting';
            app.TaucPlot.Position = [408 205 186 22];

            % Create SigmoidPlot
            app.SigmoidPlot = uicheckbox(app.UIFigure);
            app.SigmoidPlot.Text = '       Sigmoid Fit';
            app.SigmoidPlot.Position = [408 167 186 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = AutomatedApp

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end