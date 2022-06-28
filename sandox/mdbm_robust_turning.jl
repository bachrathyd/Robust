
using MDBM
using Plots


ζ=0.02::Float64
c=0.005::Float64
function D(Ω::Float64,kw::Float64,ω::Float64)::Vector{Float64}
    λ=1.0im*ω;
    τ=2.0π/Ω;
    D=λ.^2+(2*ζ+c*τ)*λ+1+(kw+1)-kw*exp(-λ*τ)
    [real(D),imag(D)]
end

#function Drs()

D(0.8,0.2,0.6)
axis=[Axis(0.25:0.2:1.5,:Ω),
Axis(-0.0:0.1:0.5,:kw),
Axis(-0.0:0.2:1.0π,:ω)]

iteration=3;
@time stab_border_points=getinterpolatedsolution(solve!(MDBM_Problem(D,axis),iteration));

scatter(stab_border_points[1],stab_border_points[2],#stab_border_points[3],
xlabel="Ω",ylabel="kw",label="",title="Stability chart of turning process with SSV",
guidefontsize=14,tickfont = font(10),markersize=2,markerstrokewidth=0)

