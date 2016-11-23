% Simple finite difference solution to the 1D wave equation. Gradients in
% time and space are both calculated using a 2nd-order accurate central
% difference scheme. The initial pressure field is set to a Gaussian.
clear all
close all

% set the literals (hard-coded numbers used in the script)
Nx = 100;       % number of grid points
Ny = 200;       % number of grid points
dx = 1e-3;      % x grid spacing (m)
dy = 1e-3;      % y grid spacing (m)
c0 = 1500;      % wave speed (m/s)
Nt = 170;       % number of time steps

% create the grid axis
x = (1:Nx)*dx;
y = (1:Ny)*dy;

% set the position of the source
x_pos = (Nx/2)*dx;
y_pos = (Ny/20)*dy;


% set the size of the time step to the stability limit, where in 1D the
% stability limit is given by dt <= dx / c0
dt = dx/(c0 * (2)^(0.5));

% set the initial xave to be a gaussian
variance = (2*dx)^2;
gaussian_x = exp( -(x - x_pos).^2 / (2 * variance));
gaussian_y = exp( -(y - y_pos).^2 / (2 * variance));
p_n = gaussian_x' * gaussian_y;

% set pressure at (n - 1) to be equal to n (this implicitly sets the
% initial particle velocity to be zero)
p_nm1 = p_n;

% preallocate the pressure at (n + 1) (this is updated during the time
% loop)
p_np1 = zeros(size(p_n));

%Create a barrier at column 30 with 2 slots of width 3dx
p_np1(1:35,90) = 0; 
p_np1(40:56,90) = 0;
p_np1(65:end,90) = 0; 


% open a new figure window
figure;

out_plot = [];

% calculate pressure in a loop
for n = 1:Nt
    
    % -----------
    % CALCULATION
    % -----------
    
    % calculate the new value for the pressure
    p_np1(2:end-1,2:end-1) = (2*p_n(2:end-1,2:end-1) - p_nm1(2:end-1,2:end-1) + ...
        (c0*dt/dx)^2 * ( p_n(1:end-2,2:end-1) - 2*p_n(2:end-1,2:end-1) + p_n(3:end,2:end-1)+ ...
        p_n(2:end-1,1:end-2) - 2*p_n(2:end-1,2:end-1) + p_n(2:end-1,3:end)))*exp(n*-1e-05);
    %Ensure the barrier stays for every iteration
    p_np1(1:35,90) = 0; 
    p_np1(40:56,90) = 0;
    p_np1(65:end,90) = 0; 
    
    %Get rid of the 'rebound' wave (based on iteration - should be a more objective way of doing this)
 if n > 140
    p_np1(:,1:89) = 0;
end
        
 
    
    % copy the value at n to (n - 1)
    p_nm1 = p_n;
    
    % copy the pressure at (n + 1) to n
    p_n = p_np1;
    
    %Plot the wave at iteration n=170 on the 'screen' at y=140

   
    
    
    % -----------
    % PLOTTING
    % -----------
    
    % plot the pressure field
    %plot(p_np1, 'k-');
    surf(p_np1);
   
    % set the limits on the z-axis
    set(gca, 'ZLim', [-0.25, 1]);
    
    % add a title with the iteration number
    title(['n = ' num2str(n)]);
    
    % force the plot to update
    drawnow;
    
    % briefly pause before continuing the loop
    pause(0.01)
        
end
figure; plot(p_np1(:,140))