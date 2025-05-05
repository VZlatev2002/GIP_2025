function cost = computeRMS(parameter_vals, parameter_names, model_name)
    load_system(model_name)
    for n = 1:height(parameter_names)
        set_param(parameter_names(n,1), parameter_names(n,2), num2str(parameter_vals(n)));
    end
    try
        simOut = sim(model_name);
    catch
        cost = 999;
        return
    end
    time_interp = linspace(0, simOut.acceleration.Time(end), 10000);
    
    acceleration_interp = interp1(simOut.acceleration.Time, simOut.acceleration.Data, time_interp);
    
    %filtered_acceleration = lsim(tf([50 500], [1, 50, 1200]), acceleration_interp, time_interp);
    
    cost = rms(acceleration_interp);
    for n = 1:height(parameter_names)
        text = join([parameter_names(n,1),": ", num2str(parameter_vals(n))]);
        disp(text)
    end
    text = join(["RMS :", cost]);
    disp(text)
    disp(' ')
    close_system(model_name, 0)
end

% function [optimal_params, optimal_rms] = NL0_optimise()
%     lower_bounds = [500, 1];
%     upper_bounds = [1000, 50];
%     initial_guess = [750, 25];
%     parameter_names = ["NL0/c_max", 'constant';
%                         "NL0/c_min", 'constant'];
%     model_name = 'NL0';
%     f = @(x)computeRMS(x,parameter_names, model_name);
%     options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
%     [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
% end
% 
% function [optimal_params, optimal_rms] = NL1_optimise()
%     lower_bounds = [1, 1, 100, 1];
%     upper_bounds = [5000, 5000, 100000, 50];
%     initial_guess = [600, 50, 10000, 50];
%     parameter_names = ["NL1/c_max", 'constant';
%                         "NL1/c_min", 'constant';
%                         "NL1/spring", 'spr_rate';
%                         "NL1/Translational Inerter", 'B'];
%     model_name = 'NL1';
%     f = @(x)computeRMS(x,parameter_names, model_name);
%     options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
%     [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
% end
%     lower_bounds = [0.79];
%     upper_bounds = [30];
%     initial_guess = [20];
%     parameter_names = ["L0/orifice1", "restriction_area"];
%     model_name = 'L0';
%     f = @(x)computeRMS(x,parameter_names, model_name);
%     options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
%     [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
% end

function [optimal_params, optimal_rms] = HNL0_optimise()
    lower_bounds = [0.79];
    upper_bounds = [30];
    initial_guess = [4];
    parameter_names = ["HNL0/orifice1", "restriction_area"];
    model_name = 'HNL0';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

function [optimal_params, optimal_rms] = HNL0b_optimise()
    lower_bounds = [0.79, 0.79, 0.79];
    upper_bounds = [30, 30, 30];
    initial_guess = [1, 10, 10];
    parameter_names = ["HNL0b/orifice1", "restriction_area";
                        "HNL0b/PRV1", "area_max";
                        "HNL0b/PRV2", "area_max"];
    model_name = 'HNL0b.slx';
    f = @(x)computeRMS(x,parameter_names, model_name);

    A = [0, 1, -1];
    b = [0];
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], A, b, lower_bounds, upper_bounds, [], options);
end

function [optimal_params, optimal_rms] = HNL1_optimise()
    lb = [1.25, 0.25, 12500];
    ub = [6, 10, 20000];
    initial_guess = [2.26, 0.41, 15000];
    % initial_guess = lower_bounds + rand(1, numel(lower_bounds)) .* (upper_bounds - lower_bounds);
    parameter_names = ["HNL1/orifice1", "restriction_area";
                        "HNL1/pipe", "length";
                        "HNL1/spring", "spr_rate"];
    model_name = 'HNL1';
    f = @(x)computeRMS(x,parameter_names, model_name);
    
    % options = optimoptions('ga', ...
    % 'UseParallel', true, ...       % Enable parallel evaluations
    % 'PopulationSize', 10, ...      % Size of the population
    % 'MaxGenerations', 100, ...     % Number of generations
    % 'Display', 'iter'); 
    % [optimal_params, optimal_rms] = ga(f, 3, [], [], [], [], lb, ub, [], options);
    % 
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lb, ub, [], options);
end

function [optimal_params, optimal_rms] = HNL1b_optimise()
    lower_bounds = [0.5, 0.1, 12500, 0.1, 0.1];
    upper_bounds = [50, 10, 17500, 30, 30];
    initial_guess = [1.5, 6, 13492, 10, 10];
    parameter_names = ["HNL1b/orifice1", "restriction_area";
                        "HNL1b/pipe", "length";
                        "HNL1b/spring", "spr_rate";
                        "HNL1b/PRV1", "area_max";
                        "HNL1b/PRV2", "area_max"];
    model_name = 'HNL1b.slx';
    f = @(x)computeRMS(x,parameter_names, model_name);

    A = [0, 0, 0, 1, -1];
    b = [0];
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], A, b, lower_bounds, upper_bounds, [], options);
end

% function [optimal_params, optimal_rms] = HNL0p_optimise()
%     lower_bounds = [0.79];
%     upper_bounds = [40];
%     initial_guess = [20];
%     parameter_names = ["HNL0p/orifice1", "restriction_area"];
%     model_name = 'HNL0p';
%     f = @(x)computeRMS(x,parameter_names, model_name);
%     options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
%     [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
% end

[HNL0_params, HNL0_J] = HNL0_optimise();
% [HNL0b_params, HNL0b_J] = HNL0b_optimise();
% [HNL1_params, HNL1_J] = HNL1_optimise();

% [HNL1b_params, HNL1b_J] = HNL1b_optimise();