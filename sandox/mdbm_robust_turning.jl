
using MDBM
using Plots
using ForwardDiff

ζ=0.02::Float64
#c=0.005::Float64
#function D(Ω::Float64,kw::Float64,ω::Float64,c::Float64)::Float64#::ComplexF64
function D(Ω,kw,ω,c)
    λ=1.0im*ω;
    τ=2.0π/Ω;
    D=λ.^2+(2*ζ+c*τ)*λ+1+(kw+1)-kw*exp(-λ*τ)
    #return D
    #vcat(real.(D), imag.(D))
    [real(D), imag(D)]
end

Dvec=xvec -> D(xvec...)



D(0.8,0.2,0.6,0.1)
Dvec([0.8,0.2,0.6,0.1])


g = x -> ForwardDiff.jacobian(Dvec, x)
g([0.8,0.2,0.6,0.1])

function Drs(x)
    g = x -> ForwardDiff.jacobian(Dvec, x)
    ComlexRoubust=imag(conj(g(x)[1,3]+1.0im*g(x)[2,3]) *(g(x)[1,:]+1.0im*g(x)[2,:]))
    vcat(Dvec(x),    ComlexRoubust[4] )
end

Drs([0.8,0.2,0.6,0.1])

#--------------------------



#------------------ Cmin --------------

axis=[Axis(0.25:0.2:1.5,:Ω),
Axis(-0.0:0.1:1.5,:kw),
Axis(-0.0:0.2:1.0π,:ω)]

iteration=2;
@time stab_border_points=getinterpolatedsolution(solve!(MDBM_Problem((x,y,z) ->D(x,y,z,0.0),axis),iteration));

scatter(stab_border_points[1],stab_border_points[2],#stab_border_points[3],
xlabel="Ω",ylabel="kw",label="",title="Stability chart of turning process with SSV",
guidefontsize=14,tickfont = font(10),markersize=2,markerstrokewidth=0)


#------------------ Cmax --------------
iteration=2;
@time stab_border_points=getinterpolatedsolution(solve!(MDBM_Problem((x,y,z) ->D(x,y,z,0.02),axis),iteration));

scatter!(stab_border_points[1],stab_border_points[2],#stab_border_points[3],
xlabel="Ω",ylabel="kw",label="",title="Stability chart of turning process with SSV",
guidefontsize=14,tickfont = font(10),markersize=2,markerstrokewidth=0)

#-----------
axis=[Axis(0.25:0.2:1.5,:Ω),
Axis(-0.0:0.1:1.5,:kw),
Axis(-0.0:0.2:1.0π,:ω),
Axis(0.00:0.01:0.02,:c)]

iteration=1;
@time stab_border_points=getinterpolatedsolution(solve!(MDBM_Problem((x,y,z,c) ->Drs([x,y,z,c]),axis),iteration));

scatter(stab_border_points[1],stab_border_points[2],#stab_border_points[3],
xlabel="Ω",ylabel="kw",label="",title="Stability chart of turning process with SSV",
guidefontsize=14,tickfont = font(10),markersize=2,markerstrokewidth=0)
