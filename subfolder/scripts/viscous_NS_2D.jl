using Plots

@views function acoustic_2D(;do_visu=false)
    # Physics
    Lx, Ly = 10.0, 10.0
    ρ      = 1.0
    μ      = 1.0
    K      = 1.0
    ttot   = 5.0
    # Numerics
    nx, ny = 128, 128
    nout   = 10
    # Derived numerics
    dx, dy = Lx/nx, Ly/ny
    dt     = min(dx,dy)^2/(K + 4/3*μ)/ρ/4.1

    nt     = cld(ttot, dt)
    xc, yc = LinRange(dx/2, Lx-dx/2, nx), LinRange(dy/2, Ly-dy/2, ny)
    # Array initialisation
    P      =  exp.(.-(xc .- Lx/2).^2 .-(yc' .- Ly/2).^2)
    dPdt   = zeros(Float64, nx,  ny)
    dVxdt  = zeros(Float64, nx-1,ny-2)
    dVydt  = zeros(Float64, nx-2,ny-1)
    Vx     = zeros(Float64, nx+1,ny)
    Vy     = zeros(Float64, nx,ny+1)
    ∇V     = zeros(Float64, ny, ny)
    τxx    = zeros(Float64, nx, ny)
    τyy    = zeros(Float64, nx, ny)
    τxy    = zeros(Float64, nx, ny)
    
    
    # Time loop
    for it = 1:nt
        dPdt .= .-K.* ∇V
        
        P .= P .+ dt.*dPdt

        τxx .= (2.0 .* μ .* (diff(Vx, dims=1)/dx .- 1.0/3.0 .* ∇V))
        τyy .= (2.0 .* μ .* (diff(Vy, dims=2)/dy .- 1.0/3.0 .* ∇V))
        τxy .= (0.5 .*(diff(Vx, dims=1)/dy + diff(Vy, dims=2)/dx))
        
        
        dVxdt.= .-1.0./ρ.*(diff(P[:,2:end-1],dims=1)./dx .- diff(τxx[:, 2:end-1], dims=1)./dx .-  diff(τxy[:, 2:end-1], dims=1)/dx)
        dVydt.= .-1.0./ρ.*(diff(P[2:end-1,:],dims=2)./dy .- diff(τyy[2:end-1, :], dims=2)./dy .-  diff(τxy[2:end-1,:], dims=2)/dy)
        Vx[2:end-1, 2:end-1]   .= Vx[2:end-1, 2:end-1] .+ dt.*dVxdt 
        Vy[2:end-1, 2:end-1]   .= Vy[2:end-1, 2:end-1] .+ dt.*dVydt
        ∇V .= diff(Vx, dims=1)./dx .+ diff(Vy, dims=2)./dy
        
        do_visu = false
        if (it % nout == 0 && do_visu)

            opts = (aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), clims=(-0.25, 0.25), c=:davos, xlabel="Lx", ylabel="Ly", title="P time = $(round(it*dt, sigdigits=3))")
            h1 = heatmap(xc, yc, P'; opts...)
            opts = (aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:davos, xlabel="Lx", ylabel="Ly", title="Vx time = $(round(it*dt, sigdigits=3))")
            h2 = heatmap(xc, yc, Vx[1:end-1, :]'; opts...)
            opts = (aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:davos, xlabel="Lx", ylabel="Ly", title="tau_xx time = $(round(it*dt, sigdigits=3))")
            h3 = heatmap(xc, yc, τxx'; opts...)
            opts = (aspect_ratio=1, xlims=(xc[1], xc[end]), ylims=(yc[1], yc[end]), c=:davos, xlabel="Lx", ylabel="Ly", title="tau_xy time = $(round(it*dt, sigdigits=3))")
            h4 = heatmap(xc, yc, τxy'; opts...)
            
            plot(h1, h2, h3, h4, layout=(2,2))
            
            savefig("Viscous_NS.png")
        end
    end
    return xc, P
end

