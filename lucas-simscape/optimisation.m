
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

function optimal_params = mL0_optimise()
    lower_bounds = [1];
    upper_bounds = [8000];
    initial_guess = [1000];
    parameter_names = ["mL0/damper", 'D'];
    model_name = 'mL0';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3,  'UseParallel',true);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

function optimal_params = mNL0_optimise()
    lower_bounds = [500, 1];
    upper_bounds = [1000, 50];
    initial_guess = [750, 25];
    parameter_names = ["mNL0/c_max", 'constant';
                        "mNL0/c_min", 'constant'];
    model_name = 'mNL0';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

function optimal_params = mL1_optimise()
    lower_bounds = [1, 5000, 1];
    upper_bounds = [1500, 15000, 250];
    initial_guess = [600, 10000, 50];
    parameter_names = ["mL1/damper", 'D';
                        "mL1/spring", 'spr_rate';
                        "mL1/inerter", 'B'];
    model_name = 'mL1';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3,  'UseParallel',true);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

function optimal_params = mNL1_optimise()
    lower_bounds = [1, 1, 5000, 1];
    upper_bounds = [1500, 400, 15000, 200];
    initial_guess = [600, 50, 10000, 50];
    parameter_names = ["mNL1/c_max", 'constant';
                        "mNL1/c_min", 'constant';
                        "mNL1/spring", 'spr_rate';
                        "mNL1/Translational Inerter", 'B'];
    model_name = 'mNL1';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

function optimal_params = L0_optimise()
    lower_bounds = [0.79];
    upper_bounds = [30];
    initial_guess = [20];
    parameter_names = ["L0/orifice1", "restriction_area"];
    model_name = 'L0';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

function optimal_params = NL0_optimise()
    lower_bounds = [0.79, 0.79];
    upper_bounds = [30, 100];
    initial_guess = [4, 80];
    parameter_names = ["NL0/orifice1", "restriction_area";
                        "NL0/orifice2", "restriction_area"];
    model_name = 'NL0';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

function optimal_params = L1_optimise()
    lower_bounds = [0.79, 0.1, 100];
    upper_bounds = [30, 10, 15000];
    initial_guess = [10.5, 3, 5000];
    parameter_names = ["L1/orifice1", "restriction_area";
                        "L1/pipe", "length";
                        "L1/spring", "spr_rate"];
    model_name = 'L1';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

function optimal_params = NL1_optimise()
    lower_bounds = [0.79, 0.1, 12000];
    upper_bounds = [3, 5, 13000];
    initial_guess = [1.5, 2, 12500];
    parameter_names = ["NL1/orifice1", "restriction_area";
                        "NL1/pipe", "length";
                        "NL1/spring", "spr_rate"];
    model_name = 'NL1';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end


function optimal_params = NL0_optimise_parasitic()
    lower_bounds = [0.79];
    upper_bounds = [40];
    initial_guess = [20];
    parameter_names = ["NL0_parasitic/orifice1", "restriction_area"];
    model_name = 'NL0_parasitic';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3, 'UseParallel',true);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

bdclose('all')

NL0_parasitic_params = NL0_optimise_parasitic();
% mNL0_params = mNL0_optimise();
% mL1_params = mL1_optimise();
%mNL1_params = mNL1_optimise();