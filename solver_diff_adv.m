% Simple reaction-diffusion-advection simulator
% allardlab.com

show_plots = true;

% model parameters
L = 5; % Length of domain, um
tmax = 100; % Total simulation time, s

D = 0.1; % Diffusion coefficient, um^2/s
v = 0.1; % Advection velocity, um/s

% numerical parameters
dz = 0.5; % Spatial step size, um

%%% ----------------------------- %%% 
% numerical timestep
dt = 0.1*min( [dz^2/D, dz/v] );

nz = ceil(L/dz)+1;

z_axis = 0:dz:L;

ntmax = round(tmax/dt);

% Diffusion matrix
D_matrix = ( -2*diag(ones(1,nz)) + diag(ones(1,nz-1),+1)+ diag(ones(1,nz-1),-1) );

% Boundary conditions for diffusion
% periodic BCs for now. 
% D_matrix(1,end) = -sum(D_matrix(:,end));
% D_matrix(end, 1) = -sum(D_matrix(:,1));

% no flux at top
D_matrix(end, end) = -1;
% zero "Dirichlet" at bottom (this is the default)

D_matrix = D_matrix*D/(dz^2);

% Advection matrix
v_down = v/dz*( -1*diag(ones(1,nz),0) + diag(ones(nz-1,1),+1) );
v_up     = v/dz*( -1*diag(ones(1,nz),0) + diag(ones(nz-1,1),-1) );

% Boundary conditions for advection
% implicitly, they are both zero ("Dirichlet") at the top and bottom

if false % inspect the diffusion and advection matrices, just for troubleshooting purposes
    figure(62); clf;
    subplot(1,3,1);
    spy(D_matrix);
    display(["Do any columns of this matrix sum to 1?" any(sum(D_matrix))])
    subplot(1,3,2);
    spy(v_up);
    display(["Do any columns of this matrix sum to 1?" any(sum(v_up))])
    subplot(1,3,3);
    spy(v_down);
    display(["Do any columns of this matrix sum to 1?" any(sum(v_down))])
end


% initial conditions
u_D = 1/L*ones(nz,1); % uniform distribution (total integral is unity)
u_v_up = 1/L*ones(nz,1); % uniform distribution (total integral is unity)
u_v_down = 1/L*ones(nz,1); % uniform distribution (total integral is unity)

nt = 0;
while nt < ntmax

    % Update the solution
    u_D = u_D + dt * ( D_matrix * u_D ); % forward-Euler. Can change this to backward-Euler for speed
    u_v_up = u_v_up + dt * ( v_up * u_v_up  );
    u_v_down = u_v_down + dt * (  v_down * u_v_down );

    % Currently, there are no reactions between the 3 states.

    if show_plots
        figure(1); clf; hold on; box on;
        plot(z_axis, u_D, 'b');
        plot(z_axis, u_v_up, 'g');
        plot(z_axis, u_v_down, 'r');

        title(['t= ' num2str(nt*dt, '%.1f') ' s']);
        xlabel('z');
        ylabel('u(z,t)');
        set(gca,'ylim', 1.2*[0,1/L], 'xlim', [0, L]);

    end

    nt = nt + 1;

end % finished time loop