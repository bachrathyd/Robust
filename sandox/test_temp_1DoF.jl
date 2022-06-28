using ForwardDiff

f(x::Vector) = exp(x[3]/10)+100*x[2]^2;
g = x -> ForwardDiff.gradient(f, x); # g = ∇f

X=[10,10,10]
g(X)

H(x::Vector)=cos(x[1])-0.3*x[1]
g = x -> ForwardDiff.gradient(H, x); # g = ∇f


using Plots

xmin=2.5
xmax=8
xi=xmin:0.01:xmax

Hi=[H([x,]) for x in xi]
gi=[g([x,])[1] for x in xi]
plot(xi,H.(xi))
plot!(xi,gi)



qi=-0.5π:0.01:0.5π
x=q->(xmax+xmin)/2+(xmax-xmin)/2*sin(q)
x.(qi)

function Hq(q::Vector)
    H([x(q[1]),])
end
Hq([3.0,])
gqq = q -> ForwardDiff.gradient(Hq, q); # g = ∇f
gqq([3,])




Hqi=[Hq([q,]) for q in qi]
gqqi=[gqq([q,])[1] for q in qi]
plot(qi,Hqi)
plot!(qi,gqqi)

