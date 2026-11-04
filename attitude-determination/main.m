%% main attitude determination script

clear;
close all

format longG
set(0, 'DefaultLegendLocation', 'eastoutside')
set(0, 'DefaultLineLineWidth', 1.4)
set(0, 'DefaultAxesFontSize', 12)
set(0, 'DefaultTextFontSize', 12)

if isempty(which('Q2Mat'))
    addpath(genpath('../SCT/SCTAcademic'));
end
if isempty(which('QProd'))
    addpath(genpath('../common'));
end

%% read data
xData = jsondecode(readlines('../results/xdata.json', 'EmptyLineRule','skip'))';
jDData = jsondecode(readlines('../results/jd.json', 'EmptyLineRule', 'skip'));
eclipseData = jsondecode(readlines('../results/eclipse.json', 'EmptyLineRule', 'skip'));

% simulation timeframe and timestep
t0 = 0;
tf = (jDData(end)-jDData(1))*86400;
dt = tf/(length(jDData)-1);
t = t0:dt:tf;
nSim = length(t);

PLOT_SENSOR_MEASUREMENTS = false;

% temperature
T = 60;

d.iR = 1:3;
d.iV = 4:6;
d.iQ = 7:10;
d.iWSat = 11:13;
d.iWRW = 14:16;

%% sensor models

% magnetometer
calM.bias = zeros(3, 1);
calM.D = zeros(3, 3);
calM.O = eye(3);
calM.sigma = (120/3)*1e-9;

mag = MagModel(calM);

% coarse sun sensors
calCSS.nSensors = 6;
calCSS.uSensors = [eye(3) -eye(3)];
calCSS.fov = deg2rad(120);
calCSS.LUT = [
    0 2.0737;
    10 2.0442;
    20 1.9592;
    30 1.8058;
    40 1.5854;
    50 1.3486;
    60 1.0999;
    ];
calCSS.LUT(:, 1) = deg2rad(calCSS.LUT(:, 1));
calCSS.sigmaTheta = deg2rad(5);
calCSS.sigmaEdge = deg2rad(10);
calCSS.sigmaDark = 1e-3;
calCSS.yLims = [0 2.4];
calCSS.T0 = 60;
calCSS.alpha = 4.31e-3;

css = CSSModel(calCSS);

% gyro
calG.bias = deg2rad([0.1; 0.1; 0.1])./3600;
calG.s = 1e-6*[1500; 1000; 1500];
calG.kU = 1e-6*[1000; 1500; 2000];
calG.kL = 1e-6*[500; 1000; 1500];
calG.sigG = repmat(sqrt(10)*1e-7, 3, 1);
calG.sigB = repmat(sqrt(10)*1e-10, 3, 1);

gyro = GyroModel(calG);

%% simulation loop

% mekf
mekfCal = MEKFGyroCal;
wPrev = xData(d.iWSat, 1);

% results
out.qList = zeros(8, nSim);
out.thetaErrList = zeros(1, nSim);
out.qErrList = zeros(4, nSim);
out.bList = zeros(6, nSim);
out.ySList = zeros(6, nSim);
out.sList = zeros(6, nSim);
out.yLitList = zeros(6, nSim);

out.xErrList = zeros(mekfCal.nStates, nSim);
out.xSigma3List = zeros(mekfCal.nStates, nSim);

tic
for i = 1:nSim
    %%% update state
    x = xData(:, i);

    r = x(d.iR);
    v = x(d.iV);
    q = x(d.iQ);
    w = x(d.iWSat);
    jD = jDData(i);

    d.jD = jD;



    %%% reference vectors (ECI)
    
    % magnetic field reference vector
    bRef = IGRFECI(r, jD); % T
    uBRef = bRef./norm(bRef);
    
    % sun reference vector
    % TODO: improve with SunV2 and compare with SunVectorECI
    uSRef = SunV1(jD, r);
    

    %%% sensor models (body frame)
    
    bMag = mag.Measurement(q, bRef, true);
    uBMag = bMag/norm(bMag);

    [ySMeas, yLit] = css.Measurement(q, uSRef, T, true);
    uSMeas = css.CalcSunVec(ySMeas, T);

    wAvg  = 0.5*(wPrev + w);
    wGyro = gyro.Measurement(wAvg, dt);
    wPrev = w;
    
    
    if ~mekfCal.initialised
        %%% TRIAD

        % TODO: investigate choice of first vector (see Wertz pg 425 and footnote)
        % TODO: investigate empirical vs analytical triad covariance

        [ATRIAD, PTRIAD] = TRIAD([uBMag, uSMeas], [uBRef, uSRef], mag.sigma/norm(bMag), 0.1);
        qEst = DCMToQ(ATRIAD);
        if (qEst(1) < 0)
            qEst = qEst.*-1;
        end

        %%% Initialise MEKF
        % initialised when 3 valid css and vectors at least 30
        % deg apart
        sinMin = 0.5;
        sinBS  = norm(cross(uBMag, uSMeas));
        if sum(yLit == 1) == 3 && sinBS > sinMin
            fprintf('MEKF initialised at %.2f min (B-sun angle %.1f deg)\n', t(i)./60, asind(sinBS))

            x0C = zeros(15, 1);
            P0C = blkdiag( ...
                PTRIAD, ...
                (0.2*pi/(3600*180))^2.*eye(3), ...
                (0.002/3)^2.*eye(3), ...
                (0.002/3)^2.*eye(3), ...
                (0.002/3)^2.*eye(3) ...
                );
            QC = blkdiag(diag([gyro.sigG.^2; gyro.sigB.^2]), zeros(9));

            mekfCal.Initialise(x0C, qEst, P0C, QC, gyro, mag, css);
        end
    else
        %%% Step MEKF
        mekfCal.Propagate(dt, wGyro);
        mekfCal.Update(q, bRef, uSRef, T);

        qEst = mekfCal.qEst;
        if (qEst(1) < 0)
            qEst = qEst.*-1;
        end
    end

    %%% Attitude error
    [thetaErr, qErr] = QAttErr(q, qEst);


    %%% save
    out.qList(1:4, i) = q;
    out.qList(5:8, i) = qEst;

    out.thetaErrList(i) = thetaErr;
    out.qErrList(:, i) = qErr;

    out.bList(1:3, i) = QToDCM(x(d.iQ))*bRef;
    out.bList(4:6, i) = bMag;

    out.ySList(:, i) = ySMeas;

    out.sList(1:3, i) = QToDCM(x(d.iQ))*uSRef;
    out.sList(4:6, i) = uSMeas;

    out.yLitList(:, i) = yLit;

    if mekfCal.initialised
        out.xErrList(:, i) = mekfCal.CalcError(q, [gyro.bias; gyro.s; gyro.kU; gyro.kL]);
        out.xSigma3List(:, i) = mekfCal.CalcSigma3();
    end


    %%% propagate
    gyro.PropagateBias(dt);
end
toc

%% plots

if PLOT_SENSOR_MEASUREMENTS
    % plot magnetic measurement
    figure('Name', 'Magnetometer');
    tl = tiledlayout(3, 1);
    tl.Title.String = "Magnetic Fields (Body)";
    tl.Title.FontWeight = "bold";
    
    out.bList = out.bList.*1e9;
    
    nexttile
    plot(t./60, out.bList(1, :), '-')
    hold on
    plot(t./60, out.bList(4, :), 'x')
    hold off
    grid on
    xticklabels({})
    ylabel('b_x (nT)')
    ylim('padded')
    
    legend('IGRF', 'Mag')
    
    nexttile
    plot(t./60, out.bList(2, :), '-')
    hold on
    plot(t./60, out.bList(5, :), 'x')
    hold off
    grid on
    xticklabels({})
    ylabel('b_y (nT)')
    ylim('padded')
    
    nexttile
    plot(t./60, out.bList(3, :), '-')
    hold on
    plot(t./60, out.bList(6, :), 'x')
    hold off
    grid on
    ylabel('b_z (nT)')
    ylim('padded')
    
    xlabel('t (min)')
    
    % plot sun measurements
    figure('Name', 'Sun Sensors');
    tl = tiledlayout(4, 1);
    tl.Title.FontWeight = "bold";
    
    nexttile
    colororder(lines(6))
    col = lines(6);
    box on
    hold on
    for i = 1:nSim
        for j = 1:css.nSensors
            litStatus = out.yLitList(j, i);
    
            if litStatus == 1 % in fov
                spec = 'o';
            elseif litStatus == 2 % outside fov
                spec = 'x';
            elseif litStatus == 0 % unlit
                spec = '.';
            end
    
            plot(t(i)./60, out.ySList(j, i), spec, 'Color', col(j, :))
        end
    end
    hold off
    title('Sun Sensor Measurements')
    grid on
    xticklabels({})
    ylabel('y_s (V)')
    ylim('padded')
    hold on
    legP{1} = plot(nan);
    legP{2} = plot(nan);
    legP{3} = plot(nan);
    legP{4} = plot(nan);
    legP{5} = plot(nan);
    legP{6} = plot(nan);
    
    legP{7} = plot(nan, 'o', 'Color', 'k');
    legP{8} = plot(nan, 'x', 'Color', 'k');
    legP{9} = plot(nan, '.', 'Color', 'k');
    hold off
    legend([legP{:}], {'+x', '+y', '+z', '-x', '-y', '-z', 'lit in FOV', 'lit outside FOV', 'unlit'})
    
    nexttile
    plot(t./60, out.sList(1, :), '-')
    title("Sun Vectors (Body)")
    hold on
    plot(t./60, out.sList(4, :), 'x')
    hold off
    grid on
    xticklabels({})
    ylabel('s_x')
    ylim('padded')
    
    legend('s_r_e_f', 's_c_s_s')
    
    nexttile
    plot(t./60, out.sList(2, :), '-')
    hold on
    plot(t./60, out.sList(5, :), 'x')
    hold off
    grid on
    xticklabels({})
    ylabel('s_y')
    ylim('padded')
    
    nexttile
    plot(t./60, out.sList(3, :), '-')
    hold on
    plot(t./60, out.sList(6, :), 'x')
    hold off
    grid on
    ylabel('s_z')
    ylim('padded')
    
    xlabel('t (min)')
end

% attitude plot
figure('Name', 'Estimated Attitude')
colororder(lines(4))

plot(t./60, out.qList(1:4, :), '-', 'DisplayName', 'q_t_r_u_e')
hold on
plot(t./60, out.qList(5:8, :), 'x','DisplayName', 'q_e_s_t')
hold off
title('Estimated Quaternion')
ylabel('q')
xlabel('t (min)')
legend()
grid on

% attitude error plots
figure('Name', 'Attitude Error')
tiledlayout(2, 1)

nexttile
plot(t./60, out.qErrList)
title('Error Quaternion')
ylabel('q error')
xlabel('t (min)')
grid on
ylim('padded')
legend('q_s', 'q_x', 'q_y', 'q_z')

nexttile
plot(t./60, rad2deg(out.thetaErrList), '-', 'MarkerSize', 8)
title('Angular Error')
ylabel('\theta error (deg)')
xlabel('t (min)')
grid on
yscale('log')

hold on
indices = find(sum(out.yLitList == 1) == 3, 1);
for i = 1:length(indices)
    xline(t(indices(i))./60, 'r--', 'LineWidth', 1.5)
end
hold off

mekfCal.Plot(t, out.xErrList, out.xSigma3List)