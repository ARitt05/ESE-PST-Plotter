classdef ESEapp < matlab.apps.AppBase
    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                  matlab.ui.Figure
        deltaLabel                matlab.ui.control.Label
        YLabel                    matlab.ui.control.Label
        EigenvaluesYXYdeltaLabel  matlab.ui.control.Label
        LabelY                    matlab.ui.control.Label
        LabelX                    matlab.ui.control.Label
        Slider                    matlab.ui.control.Slider
        SliderLabel               matlab.ui.control.Label
        Slider_2                  matlab.ui.control.Slider
        xLabel                    matlab.ui.control.Label
        MainAxes                  matlab.ui.control.UIAxes
    end
    
    properties (Access = private)
        Xvalue = 0
        Yvalue = 0
    end
    methods (Access = private)
        function CalcEngine(app)
            
            %% Define initial params
            y = round(app.Slider.Value);
            x = y + round(app.Slider_2.Value);
            app.LabelX.Text = sprintf('X = %d', x);
            app.LabelY.Text = sprintf('Y = %d', y);
            N = 5;   % Matrix dimension
            T0 = pi; % Initial PST
            
            %% Jacobi Matrix Setup
            J = Jacobi_5x5_Matrix(x,y,0);
            disp(J)
            
            %% Inner Product <e^(-iJt)e0, e0> setup
            e0 = zeros(N,1);
            e0(1) = 1; % First element of standard basis
            e4 = zeros(N,1);
            e4(5)=1;
            
            % ESE/PST Time intervals
            t = linspace(0,10,2000);
            func = zeros(size(t));
            func2 = func;
            for k = 1:length(t)
                E = expm(-1i*J*t(k));
                func(k) = e0'*E*e0;
                func2(k) = e4'*E*e0;
            end
            
            %% Plotting
            plot(app.MainAxes,t,abs(func), t, abs(func2), 'LineWidth', 0.5);
            legend(app.MainAxes,{'ESE', 'PST'})
            
            xlabel(app.MainAxes,'t')
            ylabel(app.MainAxes,'Amplitude')
            title(app.MainAxes, 'Time Evolution of Quantum States')
            % grid on;
            function J = Jacobi_5x5_Matrix(x,y,c)
                a = c*ones(5,1); % Center Diagonal
                b = (sqrt((x^2-y^2)/2)); % Center inputs of b0
                b0 = [y,b,b,y]; % Upper and Lower Diagonal
                J = diag(a)+diag(b0,1)+diag(b0,-1);
            end
            % Rules for finding points of ESE and PST
            ESE = abs(func);
            PST = abs(func2);
            
            tol = 1e-2;
            
            idx_ESE = ESE < tol;
            [~, idx_ESE] = findpeaks(-ESE);
            idx_ESE = idx_ESE(ESE(idx_ESE) < tol);
            
            m = t(idx_ESE) < T0;
            idx_ESE = idx_ESE(m);
            t_ESE = t(idx_ESE);
            [~,idx_PST] = findpeaks(PST);
            idx_PST = idx_PST(abs(PST(idx_PST) - 1) < tol);
            if ~isempty(idx_PST)
                first_PST = idx_PST(1);
                t_PST = t(first_PST);
            else
                t_PST = NaN;
            end
            
            if ~isnan(t_PST)
                m = t(idx_ESE) < t_PST;
                idx_ESE = idx_ESE(m);
                t_ESE = t(idx_ESE);
            else
                m = t(idx_ESE) < T0;
                idx_ESE = idx_ESE(m);
                t_ESE = t(idx_ESE);
            end
            
            cla(app.MainAxes)
            hold(app.MainAxes, 'on')
            
            % Displays points of ESE only when x is multiple of y
            if mod(x,y) == 0
                plot(app.MainAxes, t, ESE, 'LineWidth', 0.5, ...
                'DisplayName','$|\langle e^{-iJt} \bf{e_0} , \bf{e_0} \rangle|$')
           plot(app.MainAxes, t, PST, 'LineWidth', 0.5, ...
                'Disp','$|(\langle e^{-iJt} \bf{e_0} , \bf{e_4} \rangle|$')
                if ~isnan(t_PST)
                    yl = ylim(app.MainAxes);
                    plot(app.MainAxes, [t_PST t_PST], yl, 'r--', 'LineWidth', 1.5, 'DisplayName', 'First Instance of PST')
                end
            else
            plot(app.MainAxes, t, ESE, 'LineWidth', 0.5, ...
                'DisplayName','$|\langle e^{-iJt} \bf{e_0} , \bf{e_0} \rangle|$')
           plot(app.MainAxes, t, PST, 'LineWidth', 0.5, ...
                'Disp','$|(\langle e^{-iJt} \bf{e_0} , \bf{e_4} \rangle|$')
                if ~isnan(t_PST)
            plot(app.MainAxes, t_ESE, ESE(idx_ESE), 'bo', 'LineWidth', 1.2, 'DisplayName', 'ESE before first PST')
                end
                if ~isnan(t_PST)
                    yl = ylim(app.MainAxes);
                    plot(app.MainAxes, [t_PST t_PST], yl, 'r--', 'LineWidth', 1.5, 'DisplayName', 'First Instance of PST')
                end
            end
            legend(app.MainAxes, 'Location', 'westoutside', 'Interpreter', 'latex', 'AutoUpdate', 'on')
            end
    end
    % Callbacks that handle component events
    methods (Access = private)
        % Value changed function: Slider
        function SliderValueChanged(app, event)
            %app.Yvalue = app.Slider.Value;
            %app.YvalueLabel.Text = sprintf("Y Value: %.3f", app.Yvalue);
            CalcEngine(app)
        end
        % Value changed function: Slider_2
        function Slider_2ValueChanged(app, event)
            %app.Xvalue = app.Slider_2.Value;
            %app.XvalueLabel.Text = sprintf("X Value: %.3f", app.Xvalue);
            CalcEngine(app)
        end
        % Callback function
        function ArbitraryXYCheckBoxValueChanged(app, event)
            %value = app.ArbitraryXYCheckBox.Value;
            CalcEngine(app)
        end
    end
    % Component initialization
    methods (Access = private)
        % Create UIFigure and components
        function createComponents(app)
            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100.2 100.2 774 480];
            app.UIFigure.Name = 'MATLAB App';
            % Create MainAxes
            app.MainAxes = uiaxes(app.UIFigure);
            title(app.MainAxes, 'Title')
            xlabel(app.MainAxes, 'X')
            ylabel(app.MainAxes, 'Y')
            zlabel(app.MainAxes, 'Z')
            app.MainAxes.Position = [1 231 774 250];
            % Create xLabel
            app.xLabel = uilabel(app.UIFigure);
            app.xLabel.HorizontalAlignment = 'right';
            app.xLabel.Position = [648 39 25 22];
            app.xLabel.Text = '';
            % Create Slider_2
            app.Slider_2 = uislider(app.UIFigure);
            app.Slider_2.Limits = [1 19];
            app.Slider_2.MajorTicks = [1 3 5 7 9 11 13 15 17 19];
            app.Slider_2.Orientation = 'vertical';
            app.Slider_2.ValueChangedFcn = createCallbackFcn(app, @Slider_2ValueChanged, true);
            app.Slider_2.Step = 2;
            app.Slider_2.Position = [695 48 3 133];
            app.Slider_2.Value = 1;
            % Create SliderLabel
            app.SliderLabel = uilabel(app.UIFigure);
            app.SliderLabel.HorizontalAlignment = 'right';
            app.SliderLabel.Position = [534 42 25 22];
            app.SliderLabel.Text = '';
            % Create Slider
            app.Slider = uislider(app.UIFigure);
            app.Slider.Limits = [1 19];
            app.Slider.MajorTicks = [1 3 5 7 9 11 13 15 17 19];
            app.Slider.Orientation = 'vertical';
            app.Slider.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.Slider.Step = 2;
            app.Slider.Position = [581 51 13 127];
            app.Slider.Value = 1;
            % Create LabelX
            app.LabelX = uilabel(app.UIFigure);
            app.LabelX.Position = [25 160 77 22];
            app.LabelX.Text = 'X';
            % Create LabelY
            app.LabelY = uilabel(app.UIFigure);
            app.LabelY.Position = [24 181 77 22];
            app.LabelY.Text = 'Y';
            % Create EigenvaluesYXYdeltaLabel
            app.EigenvaluesYXYdeltaLabel = uilabel(app.UIFigure);
            app.EigenvaluesYXYdeltaLabel.Interpreter = 'latex';
            app.EigenvaluesYXYdeltaLabel.Position = [25 111 188 50];
            app.EigenvaluesYXYdeltaLabel.Text = {'Eigenvalues: '; '$Y$ & $X = Y + \delta $'};
            % Create YLabel
            app.YLabel = uilabel(app.UIFigure);
            app.YLabel.FontSize = 14;
            app.YLabel.Interpreter = 'latex';
            app.YLabel.Position = [548 48 25 22];
            app.YLabel.Text = '$Y$';
            % Create deltaLabel
            app.deltaLabel = uilabel(app.UIFigure);
            app.deltaLabel.FontSize = 14;
            app.deltaLabel.Interpreter = 'latex';
            app.deltaLabel.Position = [665 48 25 22];
            app.deltaLabel.Text = '$\delta$';
            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end
    % App creation and deletion
    methods (Access = public)
        % Construct app
        function app = ESEapp_Working
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
 
