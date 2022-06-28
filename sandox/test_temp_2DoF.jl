using ForwardDiff

H(x::Vector)=cos(x[1])-0.3*x[1]+cos(x[2])+1.3*(x[1]-5)^2+0.3*(x[2]-4)^2
g = x -> ForwardDiff.gradient(H, x); # g = ∇f


using Plots

xmin=2.5
xmax=8

ymin=2.5
ymax=7
RA=[robustaxis(xmin,xmax),robustaxis(ymin,ymax)]


xi==cyclicaxis(RA[1],20)
yi=cyclicaxis(RA[2],30)
Hi=[H([x,y,]) for y in yi, x in xi]
plot(,yi,Hi,st=:surface)
plot(xi,yi,Hi)

gi1=[g([x,y,])[1]  for y in , x in xi]
gi2=[g([x,y,])[2]  for y in yi, x in xi]
#-------------------------------

qi=-0.6π:0.01:0.6π
#qi=-0.5π:0.01:0.5π
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

