using ForwardDiff
using Plots

H(x::Vector)=cos(x[1])-0.3*x[1]+cos(x[2])+1.3*(x[1]-5)^2+0.3*(x[2]-4)^2-0.15*x[2]*x[1]
g = x -> ForwardDiff.gradient(H, x); # g = ∇f



xmin=2.5
xmax=8

ymin=2.5
ymax=6
RA=[robustaxis(xmin,xmax),robustaxis(ymin,ymax)]

xi=cyclicaxis(RA[1],120);
yi=cyclicaxis(RA[2],130);

Hi=[H([x,y,]) for y in yi, x in xi];


E=-0.01
gi1=[g([x,y,])[1] for y in yi, x in xi];
gi2=[g([x,y,])[2]  for y in yi, x in xi];
gi1=[g([x,y,])[1] *(x-xmax-E)*(x-xmin+E)*(y-ymax-E)*(y-ymin+E) for y in yi, x in xi];
gi2=[g([x,y,])[2] *(x-xmax-E)*(x-xmin+E)*(y-ymax-E)*(y-ymin+E)  for y in yi, x in xi];
gi1=[g([x,y,])[1] *(y-ymax-E)*(y-ymin+E) for y in yi, x in xi];
gi2=[g([x,y,])[2] *(x-xmax-E)*(x-xmin+E) for y in yi, x in xi];
gi1=[g([x,y,])[1] *(x-xmax-E)*(x-xmin+E) for y in yi, x in xi];
gi2=[g([x,y,])[2] *(y-ymax-E)*(y-ymin+E) for y in yi, x in xi];
ginorm=[sum(g([x,y,]).^2)  for y in yi, x in xi];
ginorm=gi1.^2 .+ gi2.^2;

plot(xi,yi,Hi,st=:surface)
plot(xi,yi,Hi)

plot(xi,yi,log.(ginorm),c=:red)
contour!(xi,yi,gi1,levels=[0,], c=:blue)
contour!(xi,yi,gi2,levels=[0,], c=:green)
#-------------------------------
function Hq(q::Vector,ras::Vector{robustaxis})
    H(axismapping(q,ras))
end

q1=-0.50π:0.01:.50π
q2=-0.50π:0.02:.50π
#qi=-0.5π:0.01:0.5π

xii=[axismapping(q,RA[1]) for q in q1]
yii=[axismapping(q,RA[2]) for q in q2]

gqq = q -> ForwardDiff.gradient(x -> Hq(x,RA), q); # g = ∇f
#Hq([3.0,3.2],RA)
#gqq([3.0,0.1])



Hqi=[Hq([qi1,qi2,],RA) for qi2 in q2, qi1 in q1];
gqi1=[gqq([qi1,qi2,])[1] for qi2 in q2, qi1 in q1];
gqi2=[gqq([qi1,qi2,])[2] for qi2 in q2, qi1 in q1];
gqinorm=gqi1.^2 .+ gqi2.^2;

plot(q1,q2,Hqi,st=:surface)
plot(q1,q2,Hqi)

plot(q1,q2,log.(gqinorm))
contour!(q1,q2,gqi1,levels=[0,], c=:blue)
contour!(q1,q2,gqi2,levels=[0,], c=:green)



plot(xii,yii,Hqi,st=:surface)
plot(xii,yii,Hqi)

plot(xii,yii,log.(gqinorm))
contour!(xii,yii,gqi1,levels=[0,], c=:blue)
contour!(xii,yii,gqi2,levels=[0,], c=:green)

