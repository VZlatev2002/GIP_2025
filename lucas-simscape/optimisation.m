
function cost = computeRMS(parameter_vals, parameter_names, model_name)

    for n = 1:height(parameter_names)
        set_param(parameter_names(n,1), parameter_names(n,2), num2str(parameter_vals(n)));
    end

    simOut = sim(model_name);
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
    
end


function L0_optimise()
    lower_bounds = [0.79];
    upper_bounds = [30];
    initial_guess = [20];
    parameter_names = ["L0/orifice1", "restriction_area"];
    model_name = 'L0';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

function L0_optimise_mechanical()
    lower_bounds = [280];
    upper_bounds = [8000];
    initial_guess = [4000];
    parameter_names = ["L0_mechanical/damper", 'D'];
    model_name = 'L0_mechanical';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

%NL0 optimisation
function NL0_optimise()
    lower_bounds = [0.79, 0.79];
    upper_bounds = [30, 100];
    initial_guess = [4, 80];
    parameter_names = ["NL0/orifice1", "restriction_area";
                        "NL0/orifice2", "restriction_area"];
    model_name = 'NL0';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

function NL0_optimise_mechanical()
    lower_bounds = [0, 0];
    upper_bounds = [8000, 8000];
    initial_guess = [4000, 30];
    parameter_names = ["NL0_mechanical/c_max", 'constant';
                        "NL0_mechanical/c_min", 'constant'];
    model_name = 'NL0_mechanical';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

function L1_optimise()
    lower_bounds = [0.79, 0.1, 100];
    upper_bounds = [30, 10, 15000];
    initial_guess = [10.5, 3, 5000];
    parameter_names = ["L1/orifice1", "restriction_area";
                        "L1/pipe", "length";
                        "L1/spring", "spr_rate"];
    model_name = 'L1';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

function NL1_optimise()
    lower_bounds = [0.79, 0.79, 0.1, 100];
    upper_bounds = [30, 50, 10, 15000];
    initial_guess = [10.5, 25, 3, 5000];
    parameter_names = ["NL1/orifice1", "restriction_area";
                        "NL1/orifice2", "restriction_area";
                        "NL1/pipe", "length";
                        "NL1/spring", "spr_rate"];
    model_name = 'NL1';
    f = @(x)computeRMS(x,parameter_names, model_name);
    options = optimoptions('patternsearch', 'Display', 'iter', 'UseCompletePoll', true, 'UseCompleteSearch', true, 'PollMethod','GPSPositiveBasis2N', 'MeshTolerance',1e-3);
    [optimal_params, optimal_rms] = patternsearch(f, initial_guess, [], [], [], [], lower_bounds, upper_bounds, [], options);
end

NL1_optimise