clear
clc
close()

% This code solves the movement equation of a 2D pendulum using three
% different integration schemes
%Legend: MD = Middle-point rule, EE = explicit Euler, EI = implicit Euler



% Numerical parameters:
T = 30;        % s
deltaT = 0.01;   % s
phi_o = 0.175/2;    % Initial angle
phi_dot_o = 0;      % Initial angular velocity
StartingConditions = [phi_o; phi_dot_o];    %Initial conditions vector
EqType = "linear";  % Equation to be solved: linear or nonlinear

% Calculation
% These arrays have the following structure [time; angle; angular speed]
% They contain all pieces of information needed to describe the system in
% minimal coordinates
EESolutions = Explicit(T, deltaT, StartingConditions, EqType);
EISolutions = Implicit(T, deltaT, StartingConditions, EqType);
MPSolutions = Midpoint(T, deltaT, StartingConditions, EqType);
ANSolutions = Analytical(T, deltaT, StartingConditions, "linear"); % Here only linear


% Comparison plot: Analytical solutions vs numerical methods:
% There are two subplots: top--> analytica; bottom -->numerical
figure(1);
h(2)=subplot(2,1,2);
box on
scatter(EESolutions(:,1), EESolutions(:,2), "blue", ".")
hold on
scatter(EISolutions(:,1), EISolutions(:,2), "black", ".")
hold on
scatter(MPSolutions(:,1), MPSolutions(:,2),"magenta", ".")
xlabel("Time [s]")
ylabel("Angle [rad]")
legend(["EE", "EI", "MP"],"Location","northwest")
hold off
box on

h(1)=subplot(2,1,1);
box on
scatter(ANSolutions(:,1), ANSolutions(:,2), "green",".")
set(h(1),'xticklabel',[]);
ylabel("Angle [rad]")
legend("ANA", "Location","northwest")

% Both together:
pos=get(h,'position');
bottom=pos{2}(2);
top=pos{1}(2)+pos{1}(4);
plotspace=top-bottom;
pos{2}(4)=plotspace/2;
pos{1}(4)=plotspace/2;
pos{1}(2)=bottom+plotspace/2;

set(h(1),'position',pos{1});
set(h(2),'position',pos{2});
linkaxes([h(1), h(2)], "y")
box on

cleanfigure;
matlab2tikz('height', '\figH', 'width', '\figW', 'filename', ['Figure1', '.tex'], 'showInfo', false, 'floatformat', '%.4g');

%%

% Constructing the energy plot:

% These arrays have the structure [time; total energy]
EEEnergy = Energy(EESolutions);
EIEnergy = Energy(EISolutions);
MPEnergy = Energy(MPSolutions);
ANEnergy = Energy(ANSolutions);

figure(2);
plot(EEEnergy(:,1), EEEnergy(:,2))
hold on
plot(EIEnergy(:,1), EIEnergy(:,2))
hold on 
plot(MPEnergy(:,1), MPEnergy(:,2))
hold on
plot(ANEnergy(:,1), ANEnergy(:,2))
hold off
legend(["EE"; "EI"; "MP"; "ANA"], "Location","northwest")
colororder(["red", "green","blue", "magenta"])
xlabel("Time [s]")
ylabel("Energy [J]")
box on

% Clearing memory:
clear EESolutions EISolutions MPSolutions ANSolutions
clear EEEnergy EIEnergy MPEnergy ANEnergy

cleanfigure;
matlab2tikz('height', '\figH', 'width', '\figW', 'filename', ['Figure2', '.tex'], 'showInfo', false, 'floatformat', '%.4g');

%%


% Convergence Analysis:
deltaTArray = [0.005, 0.001, 0.0005, 0.0001, 0.00001, 0.00005, 0.000001]; %% Time steps to be tested

MSEArray = zeros(length(deltaTArray),4);                                       % This array contains the MSE at different deltaT-values
for i=1: length(deltaTArray)
    deltaT = deltaTArray(i);     % interpolate between deltaT_min and deltaT_max
    ZustandsvektorREF = Analytical(T, deltaT, StartingConditions, "linear");     % benchmark array for all other arrays

    ZustandEE = Explicit(T, deltaT, StartingConditions, "linear");
    MSE = QAError(ZustandEE, ZustandsvektorREF);
    MSEArray(i,1) = deltaT;
    MSEArray(i,2) = MSE;
    clear ZustandEE

    ZustandEI = Implicit(T, deltaT, StartingConditions, "linear");
    MSE = QAError(ZustandEI, ZustandsvektorREF);
    MSEArray(i,1) = deltaT;
    MSEArray(i,3) = MSE;
    clear ZustandEI

    ZustandMP = Midpoint(T, deltaT, StartingConditions, "linear");
    MSE = QAError(ZustandMP, ZustandsvektorREF);
    MSEArray(i,1) = deltaT;
    MSEArray(i,4) = MSE;
    clear ZustandMP

    

end


%% Regression:
EECoeff = LinRegCoeff(log(MSEArray(:,1)), log(MSEArray(:,2)));
EICoeff = LinRegCoeff(log(MSEArray(:,1)), log(MSEArray(:,3)));
MPCoeff = LinRegCoeff(log(MSEArray(:,1)), log(MSEArray(:,4)));

EELineFit = [ones(length(MSEArray(:,1)), 1) MSEArray(:,1)]*EECoeff;
EILineFit = [ones(length(MSEArray(:,1)), 1) MSEArray(:,1)]*EICoeff;
MPLineFit = [ones(length(MSEArray(:,1)), 1) MSEArray(:,1)]*MPCoeff;


% Convergence plot: log MSE vs log deltaT
figure(3);
loglog(MSEArray(:,1), MSEArray(:,2), "r")
hold on
loglog(MSEArray(:,1), MSEArray(:,3), "g")
hold on
loglog(MSEArray(:,1), MSEArray(:,4), "b")

%{
hold on
loglog(MSEArray(:,1),EELineFit, "--")
hold on
loglog(MSEArray(:,1),EILineFit, "--")
hold on
loglog(MSEArray(:,1),MPLineFit, "--")
%}

xlabel("log(deltaT) [-]")
ylabel("log(MSE) [-]")
legend(["EE", "EI", "MP"], "Location","northwest")
colororder(["#FF0000";"#00FF00";"#0000FF"])
grid on
box on
hold off

cleanfigure;
matlab2tikz('height', '\figH', 'width', '\figW', 'filename', ['Figure3', '.tex'], 'showInfo', false, 'floatformat', '%.4g');



% Experimental convergence order is given by the slope of the line
disp("Order EE = " + num2str(EECoeff(2)))
disp("Order EI = " + num2str(EICoeff(2)))
disp("Order MP = " + num2str(MPCoeff(2)))





% Functions:



% Physical Constants
function C = Constants()
    C.g = 9.81;         % m/s^2
    C.m = 1;            % kg
    C.l = 1;            % m
    C.PhConst = C.g/C.l;
end



% Convergence Analysis


function QAError = QAError(ZustandMethod, ZustandREF)       % Calculating the time averaged mean square error
    SumSQError = 0;
    for t=1:size(ZustandREF,1)
        SumSQError = (ZustandMethod(t,2) - ZustandREF(t,2))^2 + SumSQError;
    end
    QAError = sqrt(SumSQError/size(ZustandREF,1));
end


function LinRegression = LinRegCoeff(x,y)               % For linear regression y = ax + b
    X = [ones(length(x),1) x];
    LinRegression = X\y;                                % Returned array has the structure [b; a]
end



% Solvers

function Zustandsvektor = Explicit(T, deltaT, StartingConditions, EqType)
    % initialization
    t_num = ceil(T/deltaT);         % Calculating the number of integration steps
    X = StartingConditions;         % Initialize Variable
    Zustandsvektor = zeros(t_num, 3);    % This array will store the angles and velocities for each point in time
    Zustandsvektor(1,1)= 0;
    Zustandsvektor(1,2) = X(1);
    Zustandsvektor(1,3) = X(2);
    C = Constants();
    
    % linear or non-linear solver:
    if EqType == "linear"
        Zustandsmatrix = [0 1;-C.PhConst 0]; % Initializing the Zustandsmatrix
    elseif EqType == "nonlinear"
        Zustandsmatrix = [X(2); -C.g/C.l.*sin(X(1))];
    end
    
    % Progression in time
    for t=2:t_num
        % Integration:
        if EqType =="linear"
            X_new = X + (Zustandsmatrix*X).*[deltaT; deltaT];
            X = X_new;
        elseif EqType =="nonlinear"
            X_new = X + Zustandsmatrix.*deltaT;
            X = X_new;
            Zustandsmatrix = [X(2); -C.g/C.l.*sin(X(1))];
        end

        
        % Saving the computed values:
        Zustandsvektor(t,1) = (t-1).*deltaT;        % First position: time
        Zustandsvektor(t,2) = X(1,1);             % Second position: angle
        Zustandsvektor(t,3) = X(2,1);              % Third position: angular velocity     
             
    end
    end


function Zustandsvektor = Implicit(T, deltaT, StartingConditions, EqType)

        % initialization
        t_num = ceil(T/deltaT);         % Calculating the number of integration steps
        X = StartingConditions;         % Initialize Variable
        Zustandsvektor = zeros(t_num, 3);    % This array will store the angles and velocities for each point in time
        Zustandsvektor(1,1)= 0;
        Zustandsvektor(1,2) = X(1,1);
        Zustandsvektor(1,3) = X(2,1);
        C = Constants();

        if EqType == "linear"
            Zustandsmatrix = [0 1;-C.PhConst 0]; % Initializing the Zustandsmatrix
        elseif EqType == "nonlinear"
            Zustandsmatrix = [X(1); X(2)];
        end

        % Progression in time:

        for t=2:t_num
            % integration
            if EqType == "linear"
                B = -(deltaT.*Zustandsmatrix-[1 0; 0 1]);               % see Latex-document Equation 33
                X_new = linsolve(B,X);
                X = X_new;
            elseif EqType == "nonlinear"
                fun = @(x)ImpEqSystem(x, deltaT, Zustandsmatrix);       % Optimizaiton toolbox required for this functionality
                X = fsolve(fun,Zustandsmatrix);
                Zustandsmatrix(1) = X(1);
                Zustandsmatrix(2) = X(2);
            end

            % Saving the computed values:
            Zustandsvektor(t,1) = (t-1).*deltaT;        % First position: time
            Zustandsvektor(t,2) = X(1,1);             % Second position: angle
            Zustandsvektor(t,3) = X(2,1);              % Third position: angular velocity  

        end
    
end

function F = ImpEqSystem(x, deltaT, X_alt)                  % help funciton for solving the non-linear system in the implicit case
    C = Constants();
    F(1) = x(1)-x(2).*deltaT - X_alt(1);
    F(2) = x(2) - (-C.g/C.l.*sin(x(1))).*deltaT-X_alt(2);
end



function Zustandsvektor = Midpoint(T, deltaT, StartingConditions, EqType)
     % initialization
        t_num = ceil(T/deltaT);                 % Calculating the number of integration steps
        X = StartingConditions;                 % Initialize Variable
        Zustandsvektor = zeros(t_num, 3);       % This array will store the angles and velocities for each point in time
        Zustandsvektor(1,1)= 0;                 % Time in the first column
        Zustandsvektor(1,2) = X(1);             % Angle in the second column
        Zustandsvektor(1,3) = X(2);             % Angular velocity in the third column
        C = Constants();

        if EqType == "linear"
            Zustandsmatrix = [0 1;-C.PhConst 0]; % Initializing the Zustandsmatrix
        elseif EqType == "nonlinear"
            Zustandsmatrix = [X(2)+deltaT/2.*(-C.g/C.l)*sin(X(1)); -C.g/C.l.*sin(X(1)+deltaT/2.*X(2))];
        end

        % Progression in time:
        for t=2:t_num
            % integration
            if EqType == "linear"
                X = X + deltaT.*(Zustandsmatrix*(X+deltaT/2.*(Zustandsmatrix*X)));
            elseif EqType == "nonlinear"
                X = X + Zustandsmatrix.*deltaT;
                Zustandsmatrix = [X(2)+deltaT/2.*(-C.g/C.l)*sin(X(1)); -C.g/C.l.*sin(X(1)+deltaT/2.*X(2))];
            end

            % Saving the computed values:
            Zustandsvektor(t,1) = (t-1).*deltaT;        % First position: time
            Zustandsvektor(t,2) = X(1,1);             % Second position: angle
            Zustandsvektor(t,3) = X(2,1);              % Third position: angular velocity
    
        end

end
function Zustandsvektor = Analytical(T, deltaT, StartingConditions, EqType)
    % initialization
        t_num = ceil(T/deltaT);                 % Calculating the number of integration steps
        X = StartingConditions;                 % Initialize Variable
        Zustandsvektor = zeros(t_num, 3);       % This array will store the angles and velocities for each point in time
        Zustandsvektor(1,1)= 0;
        Zustandsvektor(1,2) = X(1);
        Zustandsvektor(1,3) = X(2);
        C = Constants();
        %Zustandsmatrix = [0 1;-C.PhConst 0]; % Initializing the Zustandsmatrix

        %Progression in time:
        for t=2:t_num
            % Calculation
            if EqType == "linear"
                X(1) = StartingConditions(1).*cos(sqrt(C.PhConst).*t.*deltaT);
                X(2) = -StartingConditions(1).*sqrt(C.PhConst).*sin(sqrt(C.PhConst).*t.*deltaT);
            end

            % Analytical solution for the nonlinear case could be
            % implemented by using eliptic funcitons


            % Saving the computed values:
            Zustandsvektor(t,1) = (t-1).*deltaT;
            Zustandsvektor(t,2) = X(1);
            Zustandsvektor(t,3) = X(2);
        end

end
 
% Energy Calculation
function Energyvector = Energy(Zustandsvektor)
%Initialization
C = Constants();
I = C.m.*(C.l)^2;               % Rotational Moment of intertia
Energyvector = zeros(size(Zustandsvektor(1,1), 2));

    % Calculating total energy for each point in time
    for t=1:size(Zustandsvektor,1)
        E_pot = -C.m.*C.g.*C.l.*cos(Zustandsvektor(t,2));       % Reference point is the fixing point
        E_kin = 1/2.*I.*Zustandsvektor(t,3)^2;

        Energyvector(t,1) = Zustandsvektor(t,1);    % First column: time
        Energyvector(t,2) = E_pot + E_kin;          % Second column: total energy
    end
end

