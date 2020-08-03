
using Pkg;
pkg"registry add https://github.com/FrameFunVC/FrameFunRegistry"

using Pkg;Pkg.activate(pwd())
Pkg.instantiate()
# apath=normpath(joinpath(pwd(),"."));
width = ".5\\textwidth"
height = ".25\\textwidth"

using BasisFunctions, DomainSets, PGFPlotsX, FrameFun, LaTeXStrings, DocumentPGFPlots

N = 31
f = identity
P = ExtensionFramePlatform(platform(Fourier(10)→(-2.0..2.0)),-1.0..1.0)
F = Expansion(basis(dictionary(P,N)),coefficients(Fun(f, P, N;verbose=true)))
t = PeriodicEquispacedGrid(2000,-4,4)
tt = PeriodicEquispacedGrid(1000,-2,2)
P1Dexampleframe = @pgf Axis(
    {width=width,height=height,ymin=-2.5,ymax=2.5},
    PlotInc({mark="none", style="dashed,blue"},Table([t, circshift(repeat(real.((F)(tt)),2),500)])),
    PlotInc({mark="none",style="blue,thick"},Table([tt, real.((F)(tt))])),
    PlotInc({mark="none",style="red, thick"},Table([t, t]))
    )

t = PeriodicEquispacedGrid(2000,-4,4)
tt = PeriodicEquispacedGrid(500,-1,1)
P = platform(Fourier(10,-1,1))
F = Fun(f, P, N;verbose=true)
P1Dexamplebasis = @pgf Axis(
    {width=width,height=height,ymin=-2.5,ymax=2.5},
    PlotInc({mark="none", style="dashed,blue"},Table([t, circshift(repeat(real.((F)(tt)),4),250)])),
    PlotInc({mark="none",style="blue,thick"},Table([tt, real.((F)(tt))])),
    PlotInc({mark="none",style="red, thick"},Table([t, t]))
    )

# DocumentPGFPlots.savefigs(joinpath(apath,"img","1Dexample"),P1Dexampleframe)
# DocumentPGFPlots.savefigs(joinpath(apath,"img","1Dexample_Gibbs"),P1Dexamplebasis)

using BasisFunctions, DomainSets, PGFPlotsX, FrameFun, LaTeXStrings, DocumentPGFPlots

f = exp
N1 = 20
N2 = 40
Pbasis = platform(ChebyshevT(1,-1,1))
P = ExtensionFramePlatform(platform(ChebyshevT(1,-2.,2.)), -1.0..1.0)

@time Fbasis = Fun(f, Pbasis, N2;verbose=true)
@time Fbasissmall = Fun(f, Pbasis, N1;verbose=true)

Fframe = Fun(f, P, N2;verbose=true)
Fframesmall = Fun(f, P, N1;verbose=true)

Pbasis = @pgf Axis({width=width,height=height,ymode="log",ymin=1e-18,ymax=4,ylabel=L"$|c_k|$",xlabel=L"$k$"},
    PlotInc({},Table([1:N2,abs.(coefficients(Fbasis))])),
    PlotInc({},Table([1:N1, abs.(coefficients(Fbasissmall))]))
    )

Pframe = @pgf Axis({width=width,height=height,ymode="log",ymin=1e-18,ymax=4,ylabel=L"$|c_k|$",xlabel=L"$k$"},
    PlotInc({},Table([1:N2,abs.(coefficients(Fframe))])),
    PlotInc({},Table([1:N1,abs.(coefficients(Fframesmall))]))
    )

# DocumentPGFPlots.savefigs(joinpath(apath,"img","coefsinbasis"),Pbasis)
# DocumentPGFPlots.savefigs(joinpath(apath,"img","coefsinframe"),Pframe)

f = exp
err1 = zeros(N2)
err2 = similar(err1)
for (i,N) in enumerate(1:N2)
    @show i
    Pbasis = platform(ChebyshevT(1,-1,1))
    P = ExtensionFramePlatform(platform(ChebyshevT(1,-2.,2.)), -1.0..1.0)
    Fbasis = Fun(f, Pbasis, N;threshold=1e-14,normalizedsampling=true)
    Fframe = Fun(f, P, N;threshold=1e-14,normalizedsampling=true)
    err1[i] = maxerror(f, Fbasis)
    err2[i] = maxerror(f, Fframe)
end


Pbasis = @pgf Axis({width=width,height=height,ymode="log",ymin=1e-18,ymax=4,ylabel=L"$\|f-f_N\|_{\infty}$",xlabel=L"$N$"},
    PlotInc({},Table([1:N2,err1]))
    )

Pframe = @pgf Axis({width=width,height=height,ymode="log",ymin=1e-18,ymax=4,ylabel=L"$\|f-f_N\|_{\infty}$",xlabel=L"$N$"},
    PlotInc({},Table([1:N2,err2]))
    )

# DocumentPGFPlots.savefigs(joinpath(apath,"img","errinbasis"),Pbasis)
# DocumentPGFPlots.savefigs(joinpath(apath,"img","errinframe"),Pframe)

using BasisFunctions, DomainSets, PGFPlotsX, FrameFun, LaTeXStrings, DocumentPGFPlots

width = ".5\\textwidth"
height = ".3\\textwidth"

P = ExtensionFramePlatform(platform(ChebyshevT(1)→(-2.0..2.)), -1.0..1.0)
f = exp

N  = 61
fscale = (d,i) -> 10.0^-4+abs(i)+abs(i)^2+abs(i)^3
Frough_coef = coefficients(Fun(f, P, N;threshold=1e-14,normalizedsampling=true,solverstyle=AZStyle()))
Fsmooth_coef = coefficients(Fun(f, P, N;threshold=1e-14,normalizedsampling=true,
        solverstyle=AZSmoothStyle(),scaling=fscale))
Fwa_coef = Fun(f,P,threshold=1e-14,normalizedsampling=true,δ=0,p0=4,
    adaptivestyle=:greedy,maxlength=N-1,weightedAZ=true,maxiterations=N)
Frough = Expansion(basis(dictionary(P,N)), Frough_coef)
Fsmooth = Expansion(basis(dictionary(P,N)), Fsmooth_coef)
Fwa = Expansion(basis(dictionary(P,N)), Fwa_coef)

N2 = 101
Frough_coef2 = coefficients(Fun(f, P, N2;threshold=1e-14,normalizedsampling=true,solverstyle=AZStyle()))
Fsmooth_coef2 = coefficients(Fun(f, P, N2;threshold=1e-14,normalizedsampling=true,
        solverstyle=AZSmoothStyle(),scaling=fscale))
Fwa_coef2 = Fun(f,P,threshold=1e-14,normalizedsampling=true,δ=0,p0=4,
    adaptivestyle=:greedy,maxlength=N2-1,weightedAZ=true,maxiterations=N2)
Frough2 = Expansion(basis(dictionary(P,N2)), Frough_coef2)
Fsmooth2 = Expansion(basis(dictionary(P,N2)), Fsmooth_coef2)
Fwa2 = Expansion(basis(dictionary(P,N2)), Fwa_coef2)

s1 = @pgf {style="red",mark="*",mark_options="red"}
s2 = @pgf {style="blue",mark="square*",mark_options="blue"}
s3 = @pgf {style="green",mark="diamond*",mark_options="green"}
s4 = @pgf {style="brown",mark="diamond"}

s1d = @pgf {style="red,dashed",mark="o"}
s2d = @pgf {style="blue,dashed",mark="square"}
s3d = @pgf {style="green,dashed",mark="diamond"}
s4d = @pgf {style="brown,dashed",mark="diamond"}

t = LinRange(-2,2,1000)
P1 = @pgf Axis(
    {xlabel = L"$t$",
    legend_cell_align="left",
    width=width,height=height,ymin=-5,ymax=5,legend_pos="south west"},
    PlotInc({mark="none",style="very thick",color="red"},Table([t, real((Frough).(t))])),
    PlotInc({mark="none",style="very thick",color="blue"},Table([t, real.((Fsmooth).(t))])),
    PlotInc({mark="none",style="very thick",color="green"},Table([t, real.((Fwa).(t))])),
    )

P2 = @pgf Axis(
    {xlabel=L"$t$",
    width=width,height=height,ymode="log"},
    PlotInc({mark="none",style="very thick",color="red"},Table([t, eps() .+ abs.(Frough.(t)-f.(t))])),
    PlotInc({mark="none",style="very thick",color="blue"},Table([t, eps() .+ abs.(Fsmooth.(t)-f.(t))])),
    PlotInc({mark="none",style="very thick",color="green"},Table([t, eps() .+ abs.(Fwa.(t)-f.(t))])),
)

r =sortperm(map(x->convert(Int,x),ordering(dictionary(P,N))));
r2 =sortperm(map(x->convert(Int,x),ordering(dictionary(P,N2))));

P3 = @pgf Axis(
    {xlabel=L"$k$",
    width=width,height=height,ymode="log",
        legend_pos="south west",legend_style={legend_columns=-1,fill_opacity=0.6,
    text_opacity=1},legend_cell_align={left}},
    PlotInc(s1, Table([1:N,
                abs.(Frough_coef)])),
    PlotInc(s2, Table([1:N,
                abs.(Fsmooth_coef)])),
    PlotInc(s3, Table([1:N,
                abs.(Fwa_coef)])),

    PlotInc(s1d, Table([1:N2,
                abs.(Frough_coef2)])),
    PlotInc(s2d, Table([1:N2,
                abs.(Fsmooth_coef2)])),
    PlotInc(s3d, Table([1:N2,
                abs.(Fwa_coef2)])),
)

# DocumentPGFPlots.savefigs(joinpath(apath,"img","approximation"),P1)
# DocumentPGFPlots.savefigs(joinpath(apath,"img","error"),P2)
# DocumentPGFPlots.savefigs(joinpath(apath,"img","coefficients"),P3)



using BasisFunctions, DomainSets, PGFPlotsX, FrameFun, LaTeXStrings, DocumentPGFPlots
width = ".5\\textwidth"
height = ".3\\textwidth"

f = x->1.
m_mon = 30
err_mon = zeros(m_mon)
csize_mon = zeros(m_mon)
tol = 1e-14; cutoff = 1e-5
for N in 1:m_mon
    P = ExtensionFramePlatform(platform(Fourier(1,-2,2)),-1.0..1.0)
    A = Fun(f, P, N;samplingstyle=GramStyle(),REG=regularized_SVD_solver, threshold=cutoff,
        discrete=false, rtol=tol, atol=tol)
    err_mon[N] = L2error(f,A,  atol=tol, rtol=tol)
    csize_mon[N] = norm(coefficients(A))
end

upper_mon = [NaN, (err_mon[2:end].+sqrt(cutoff)*csize_mon[2:end])...];

Pnonmonotone = @pgf Axis(
    {width=width,height=height,ymode="log",xlabel=L"$N$",
    cycle_list_name="mark list*",
    ymin=1e-5},
    PlotInc({mark_size="1pt"},Table([1:m_mon,upper_mon])),
    PlotInc({mark_size="1pt"},Table([1:m_mon,err_mon])),
)

# DocumentPGFPlots.savefigs(joinpath(apath,"img","notmonotone"),Pnonmonotone)

f = x->exp(cos(5x))
m = 30
err = zeros(m)
csize = zeros(m)
tol = 1e-14; cutoff = 1e-5
for N in 1:m
    P = ExtensionFramePlatform(platform(Fourier(1,-2,2)),-1.0..1.0)
    A = Fun(f, P, N;samplingstyle=GramStyle(),REG=regularized_SVD_solver,
        threshold=cutoff, discrete=false, rtol=tol, atol=tol)
    err[N] = L2error(f,A,  atol=tol, rtol=tol)
    csize[N] = norm(coefficients(A))
end

upper = [NaN, (err[2:end].+sqrt(cutoff)*csize[2:end])...];

P = @pgf Axis(
    {width=width,height=height,ymode="log",xlabel=L"$N$",
    cycle_list_name="mark list*",
    ymin=1e-5},
    PlotInc({mark_size="1pt"},Table([1:m,upper])),
    PlotInc({mark_size="1pt"},Table([1:m,err])),
)

# DocumentPGFPlots.savefigs(joinpath(apath,"img","decreaseofL2norm"),P)

f = x->1.
m_mon = 40
err_mon = zeros(m_mon)
csize_mon = zeros(m_mon)
tol = 1e-14; cutoff = 1e-5
for N in 1:m_mon
    P = ExtensionFramePlatform(platform(Fourier(1,-2,2)),-1.0..1.0)
    A = Fun(f, P, N;REG=regularized_SVD_solver, threshold=cutoff)
    err_mon[N] = L2error(f,A,  atol=tol, rtol=tol)
    csize_mon[N] = norm(coefficients(A))
end

upper_mon =  [NaN, (err_mon[2:end].+(cutoff)*csize_mon[2:end])...];

Pnonmonotone = @pgf Axis(
    {width=width,height=height,ymode="log",xlabel=L"$N$",
    cycle_list_name="mark list*",
    ymax=1,ymin=1e-12},
    PlotInc({mark_size="1pt"},Table([1:m_mon,upper_mon])),
    PlotInc({mark_size="1pt"},Table([1:m_mon,err_mon])),
)

# DocumentPGFPlots.savefigs(joinpath(apath,"img","notmonotone2"),Pnonmonotone)

f = x->exp(cos(5x))
m = 40
err = zeros(m)
csize = zeros(m)
tol = 1e-14; cutoff = 1e-5
for N in 1:m
    P = ExtensionFramePlatform(platform(Fourier(1,-2,2)),-1.0..1.0)
    A = Fun(f, P, N;REG=regularized_SVD_solver, threshold=cutoff)
    err[N] = L2error(f,A,  atol=tol, rtol=tol)
    csize[N] = norm(coefficients(A))
end

upper = [NaN, (err[2:end].+(cutoff)*csize[2:end])...];

P = @pgf Axis(
    {width=width,height=height,ymode="log",xlabel=L"$N$",
    cycle_list_name="mark list*",
    ymin=1e-12},
    PlotInc({mark_size="1pt"},Table([1:m,upper])),
    PlotInc({mark_size="1pt"},Table([1:m,err])),
)

# DocumentPGFPlots.savefigs(joinpath(apath,"img","decreaseofL2norm2"),P)

using FrameFun, DomainSets, LaTeXStrings,PGFPlotsX, Statistics, ProgressBars, QuadGK
width = ".5\\textwidth"
height = ".3\\textwidth"
P = ExtensionFramePlatform(platform(Fourier(1,-2,2)),-1.0..1.0)
function getmedianindex(N)
    mindex = Array{Int}(undef,size(N)[1:end-1])
    for i in CartesianIndices(size(N)[1:end-1])
        mindex[i] = findfirst(median(N[i.I...,:]) .==N[i.I...,:])
    end
    mindex
end
function getmedian(N,mindex)
    NN = Array{eltype(N)}(undef,size(mindex))
    for i in CartesianIndices(size(NN))
        NN[i] = N[i.I...,mindex[i]]
    end
    NN
end

ps = LinRange(10,500,20)
noexp = 7
Nmax = 2maximum(ps)

N1sfull = zeros(Int, length(ps), noexp)
totaltimefull = zeros(Float64, length(ps), noexp)
N2sfull = copy(N1sfull)
onetimefull = copy(totaltimefull)
ix = 1:length(ps)

for i in ProgressBar(ix)
    for n in 1:noexp
        p = ps[i]
    #         @show p
        f = x->cos(p*x)
        _,t1,_ = @timed( F1 = Fun(f, P; numrandompts = 3, adaptivestyle = OptimalStyle(),
                criterion = FNAStyle(),δ=1e-10,threshold=1e-12,maxlength = Nmax,maxiterations = 10000))
        totaltimefull[i,n] = t1
        N1sfull[i,n] = length(coefficients(F1))
        _,t2,_ = @timed(Fun(f, P,length(coefficients(F1)),threshold=1e-12))
        onetimefull[i,n] = t2
    end
end

mindex = getmedianindex(N1sfull)
N1s = getmedian(N1sfull,mindex)
N2s = getmedian(N2sfull,mindex)
totaltime = getmedian(totaltimefull,mindex)
onetime = getmedian(onetimefull,mindex);

P1 = @pgf Axis({xlabel=L"p",width=width,height=height,cycle_list_name="mark list*",},
    Plot(Table([ps[2:end],totaltime[:,1,1][2:end]./onetime[:,1,1][2:end]]))
    )

for i in ProgressBar(ix)
    for n in 1:noexp
        p = ps[i]
#         @show p
        f = x->cos(p*x)
        F2 = Fun(f, P; adaptivestyle = GreedyStyle(), numrandompts = 3,
            p0=max(N1s[i]-10,4), criterion = FNAStyle(),maxlength = Nmax,maxiterations = 10000)
        N2sfull[i,n] = length(coefficients(F2))
        @assert N2sfull[i,n] > max(N1s[i]-10,4)
    end
end
N2s = getmedian(N2sfull,mindex);

P2 = @pgf Axis({cycle_list_name="mark list*",xlabel=L"p",width=width,height=height,ylabel=L"\Delta N"},
    PlotInc(Table([ps,N1s.-N2s])),
)

# using DocumentPGFPlots
# DocumentPGFPlots.savefigs(joinpath(apath,"img","adaptiveVsOptimal"),P1)
# DocumentPGFPlots.savefigs(joinpath(apath,"img","adaptiveOptimalN"),P2)

using FrameFun, DomainSets, LaTeXStrings,PGFPlotsX, Statistics
p = IncrementalCartesianParameterPath{2}()
P = ExtensionFramePlatform(WeightedSumPlatform(platform(ChebyshevT(10,-1,1)^2), (x,y)->1.,
            (x,y)->sqrt(x^2+y^2)),.9*UnitDisk())
path = HierarchyPath(p,ProductPath(p,p))
paramP = parametrizedplatform(P, path)
function getmedianindex(N)
    mindex = Array{Int}(undef,size(N)[1:end-1])
    for i in CartesianIndices(size(N)[1:end-1])
        mindex[i] = findfirst(median(N[i.I...,:]) .==N[i.I...,:])
    end
    mindex
end
function getmedian(N,mindex)
    NN = Array{eltype(N)}(undef,size(mindex))
    for i in CartesianIndices(size(NN))
        NN[i] = N[i.I...,mindex[i]]
    end
    NN
end

ps = LinRange(0,3,10)
noexp = 1
Nmax = 3000
N1sfull = zeros(Int, length(ps), noexp)
totaltimefull = zeros(Float64, length(ps), noexp)
onetimefull = copy(totaltimefull)

for (i,p) in ProgressBar(enumerate(ps))
    for n in 1:noexp
        f = (x,y) -> cos(p*pi*(x+y)) + sqrt(x^2+y^2)*sin(1+p*pi*(x+y))
        _,t1,_ = @timed( F1 = Fun(f, paramP; numrandompts = 3,adaptivestyle = OptimalStyle(),
                criterion = FNAStyle(),δ=1e-6, threshold=1e-8, maxlength = Nmax,maxiterations = 10000))
        totaltimefull[i,n] = t1
        N1sfull[i,n] = length(coefficients(F1))
        _,t2,_ = @timed(F=Fun(f, P,tuple(dimensions(dictionary(F1))...), threshold=1e-8, ))
        onetimefull[i,n] = t2
    end
end

mindex = getmedianindex(N1sfull)
N1s = getmedian(N1sfull,mindex)
totaltime = getmedian(totaltimefull,mindex)
onetime = getmedian(onetimefull,mindex);

P1 = @pgf GroupPlot({group_style = {group_size="2 by 1",},cycle_list_name="mark list*"},
    {width=width,height=height,xlabel=L"p",ylabel=L"N",ymin=0},
    Plot(Table([ps,N1s])),
    {legend_cell_align="left",width=width,height=height,ymin=0,ymax=20,xlabel=L"p"},
    Plot(Table([ps[2:end],totaltime[2:end]./onetime[2:end]])),
)

# using DocumentPGFPlots
# DocumentPGFPlots.savefigs(joinpath(apath,"img","adaptiveVsOptimal2dwFE"),P1)

f = (x,y) -> cos(last(ps)*pi*(x+y)) + sqrt(x^2+y^2)*sin(1+last(ps)*pi*(x+y))
F = Fun(f, paramP; numrandompts = 3,adaptivestyle = OptimalStyle(),
    criterion = FNAStyle(),δ=1e-6, threshold=1e-8,maxlength = Nmax,maxiterations = 10000,verbose=false)

using Plots
x = EquispacedGrid(100,-1,1)
Plots.plot(F;c=:RdBu,size=(2*300,2*140),layout=2)
wcP = Plots.heatmap!(log10.(eps().+abs.(F(x^2)-[f(xi, yi) for xi in x, yi in x]));subplot=2,aspect_ratio=1,ticks=false)

# savefig(wcP, joinpath(apath,"img","weightedcircle"))





using BasisFunctions, FrameFun, DomainSets, PGFPlotsX, Plots, LaTeXStrings, DocumentPGFPlots

width = ".5\\textwidth"
height = ".3\\textwidth"

ymax = 1e6; ymin = 1e-15
epss = exp10.(LinRange(-12,-3,4))
T = 2; sampling_factor=2;
t = 1;

f_1 = x->exp(cos(8pi*x))
# f_2 = x->1/(x-1.01)

ns_1 = 2collect(1:10:200) .+ 1
csize_1 = zeros(Float64,length(ns_1),length(epss))
L2err_1 = zeros(Float64,length(ns_1),length(epss))
l2err_1 = zeros(Float64,length(ns_1),length(epss))
fsize_1 = zeros(Float64,length(ns_1),length(epss))
minN_1 = zeros(Int,length(epss))
P = ExtensionFramePlatform(platform(Fourier(1, -T,T)),Interval(-t,t))
for (i,n) in enumerate(ns_1)
    for (j,epsilon) in enumerate(epss)
        F, A, b, c, S, _  = approximate(f_1,P,n;solverstyle=AZStyle(),L=sampling_factor*n,threshold=epsilon)
        csize_1[i,j] = norm(c)
        L2err_1[i,j] = L2error(f_1, F; atol=1e-5,rtol=1e-5)
        b = S*f_1
        l2err_1[i,j] = norm(b-A*c)
        fsize_1[i,j] = sqrt(step(supergrid(supergrid(FrameFun.grid(element(S,1))))).*norm(b)^2)
        (minN_1[j] == 0) && (l2err_1[i,j] < epsilon*fsize_1[i,j]) && (minN_1[j] = n)
    end
end

@pgf p_1 = GroupPlot(
    {width=width,height=height,ymode="log",
    cycle_list_name="mark list*",
    xlabel=L"$N$", legend_pos="south west",legend_cell_align="left",
        ymin=1e-13,ymax=1e12,group_style={group_size={"2 by 2"},
        x_descriptions_at="edge bottom",
        y_descriptions_at="edge left",
        horizontal_sep="1em",vertical_sep="1em",

    }},
    )
@pgf for (j,epsilon) in enumerate(epss)
    push!(p_1,{})
    push!(p_1,PlotInc({},Table(ns_1,L2err_1[:,j])))
    push!(p_1,PlotInc({},Table(ns_1,l2err_1[:,j])))
    push!(p_1,PlotInc({},Table(ns_1,csize_1[:,j])))
    push!(p_1,PlotInc({style="black",mark="none"},Table(ns_1,fsize_1[:,j]))) 
    push!(p_1,PlotInc({style="black,dashed",mark="none"},Table(ns_1, epsilon*ones(size(ns_1)))))
end
p_1

# DocumentPGFPlots.savefigs(joinpath(apath,relpath,"groupplotexpcos8pix"),p_1)



P = ExtensionFramePlatform(FourierPlatform(),0.0..0.5)

f = x->exp(x) + 1e-3*cos(2pi*1000x)
F = Fun(f,P;threshold=1e-10/100,δ=1e-10,criterion=FNAStyle(),verbose=false)

δs = 10. .^ LinRange(-10,-3,20)
ϵs = 10. .^ LinRange(-10,-3,20)
ls = zeros(length(δs),length(ϵs))
ns = zeros(length(δs),length(ϵs))
es = zeros(length(δs),length(ϵs))
for (i,δ) in enumerate(δs)
    for (j,ϵ) in enumerate(ϵs)
        f = x->exp(x) + ϵ*cos(2pi*1000x)
        F = Fun(f,P;threshold=δ/100,δ=δ,criterion=FNAStyle(),verbose=false)
        N = length(coefficients(F))
        _,A,b,c,S,L = approximate(f,P,N;threshold=δ/100,δ=δ)
        ls[i,j] = length(F)
        ns[i,j] = norm(coefficients(F))
        es[i,j] = norm(A*c-b)
        @show δ, ϵ, length(F), norm(coefficients(F)), norm(A*c-b)
    end
end

opts = @pgf {
        view = (0, 90),
        xlabel=L"\log_{10}(\delta)",
        ylabel=L"\log_{10}(\sigma)",
        colorbar,
        "colormap/jet",
        width = ".25\\textwidth",
        scale_only_axis=true,
    }
p3opts = @pgf {surf,shader = "flat"}

P = @pgf GroupPlot({opts...,group_style={
        x_descriptions_at="edge bottom",
        y_descriptions_at="edge left",
        horizontal_sep="5em",group_size="2 by 1",}},
    {},Plot3(p3opts,Table(log10.(ϵs), log10.(δs),log10.(es))),
    {},Plot3(p3opts,Table(log10.(ϵs), log10.(δs),log10.(ls))),
)

# DocumentPGFPlots.savefigs(joinpath(apath,"img","svsd"), P)

f = x-> 1e6exp(cos(8pi*x))
P = ExtensionFramePlatform(FourierPlatform(),0.0..0.5)

δs = 10. .^ LinRange(-10,-3,20)
ϵs = 10. .^ LinRange(-10,-3,20)
ls = zeros(length(δs),length(ϵs))
ns = zeros(length(δs),length(ϵs))
for (i,δ) in enumerate(δs)
    for (j,ϵ) in enumerate(ϵs)
        F = Fun(f,P;threshold=ϵ,δ=δ,criterion=FNAStyle(),verbose=false)
        @show δ, ϵ, length(F), norm(coefficients(F))
        ls[i,j] = length(F)
        ns[i,j] = norm(coefficients(F))
    end
end

f2 = x-> exp(cos(8pi*x))
P2 = ExtensionFramePlatform(FourierPlatform(),0.0..0.5)

δs2 = 10. .^ LinRange(-10,-3,20)
ϵs2 = 10. .^ LinRange(-10,-3,20)
ls2 = zeros(length(δs2),length(ϵs2))
ns2 = zeros(length(δs2),length(ϵs2))
for (i,δ) in enumerate(δs2)
    for (j,ϵ) in enumerate(ϵs2)
        F = Fun(f2,P2;threshold=ϵ,δ=δ,criterion=FNAStyle(),verbose=false)
        @show δ, ϵ, length(F), norm(coefficients(F))
        ls2[i,j] = length(F)
        ns2[i,j] = norm(coefficients(F))
    end
end

opts = @pgf {
        view = (0, 90),
        xlabel=L"\log_{10}(\delta)",
        ylabel=L"\log_{10}(\epsilon)",
        colorbar,
        "colormap/jet",
        width = ".25\\textwidth",
        scale_only_axis=true,
    }
p3opts = @pgf {surf,shader = "flat"}

P = @pgf GroupPlot({opts...,group_style={
        x_descriptions_at="edge bottom",
        y_descriptions_at="edge left",
        horizontal_sep="5em",group_size="4 by 1",}},
    {L"$f_1, N_{\mathrm{opt}}$"},Plot3(p3opts,Table(log10.(ϵs), log10.(δs),log10.(ls2))),
    {L"$f_1, \|\mathbf x\|$"},Plot3(p3opts,Table(log10.(ϵs), log10.(δs),log10.(ns2))),
    {L"$f_2, N_{\mathrm{opt}}$"},Plot3(p3opts,Table(log10.(ϵs), log10.(δs),log10.(ls))),
    {L"$f_2, \|\mathbf x\|$"},Plot3(p3opts,Table(log10.(ϵs), log10.(δs),log10.(ns))),
)

# DocumentPGFPlots.savefigs(joinpath(apath,"img","evsd"), P)
