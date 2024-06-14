module solver_diff_adv

using Plots
using LinearAlgebra
using Printf
using ProgressBars

function main()
    # show_plots = true

    L = 5
    tmax = 100

    D = 0.1
    v = 0.1

    dz = 0.5

    dt = 0.1 * min(dz^2/D, dz/v)
    nz = Integer(ceil(L/dz) + 1)
    z_axis = 0:dz:L
    ntmax = round(tmax/dt)

    D_matrix = diagm(0 => repeat([-2], nz),
                     1 => repeat([1], nz-1),
                     -1 => repeat([1], nz-1)
                    )

    # no flux at top: Dirichlet?
    D_matrix[end, end] = -1

    D_matrix = D_matrix*D/dz^2

    # advection matrix
    v_down = (v/dz)*diagm(
        0 => repeat([-1], nz),
        1 => repeat([1], nz-1)
    )
    v_up = (v/dz)*diagm(
        0 => repeat([-1], nz),
        -1 => repeat([1], nz-1)
    )

    # initial conditions
    u_D = 1/L*repeat([1], nz)
    u_v_up = 1/L*repeat([1], nz)
    u_v_down = 1/L*repeat([1], nz)

    nt = 0
    ani = Animation()
    maxY = 0
    for nt = ProgressBar(1:ntmax)
        u_D = u_D + dt * ( D_matrix * u_D  )
        u_v_up = u_v_up + dt * ( v_up * u_v_up  )
        u_v_down = u_v_down + dt * ( v_down * u_v_down  )

        maxY = max(maxY, maximum(vcat(u_D, u_v_up, u_v_down)))

        plt = plot(z_axis, u_D, label="u_D")
        plot!(z_axis, u_v_up, label="u_v_up")
        plot!(z_axis, u_v_down, label="u_v_down")
        plot!(title=@sprintf("t = %.2f", nt*dt), ylim=[0,maxY])
        frame(ani, plt)

    end
    display(gif(ani, "out.mp4"))

end

main()

end
