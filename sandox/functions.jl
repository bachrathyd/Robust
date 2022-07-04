struct robustaxis
    qmin::Float64
    qmax::Float64
    Δq::Float64
    qmean::Float64
end

robustaxis(a,b)=robustaxis(a,b,(b-a)/2,(b+a)/2)

function axismapping(qi::Real,ra::robustaxis)
    ra.qmean+ra.Δq*sin(qi);
end
   
function cyclicaxis(ra::robustaxis,N::Int)
    xi=LinRange(ra.qmin,ra.qmax,N)
end


function axismapping(q::Vector,ras::Vector{robustaxis})
    [axismapping(qi,ra) for (qi,ra) in zip(q,ras)]
end