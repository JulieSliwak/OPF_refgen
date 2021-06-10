"""
    struct MatpowerLine_π_Smax <: AbstractLinkLabel

Structure descendant of AbstractLinkLabel

# Fields
- `baseMVA::Float64`
- `resistance::Float64`: coefficient rs
- `reactance::Float64`: coefficient xs
- `susceptance::Float64`: coefficient bc
- `transfo_ratio::Float64`: ratio τ
- `transfo_phase::Float64`: phase θ

"""
struct MatpowerLine_π_Smax <: AbstractLinkLabel
  baseMVA::Float64
  resistance::Float64 #rs
  reactance::Float64 #xs
  susceptance::Float64 #bc
  transfo_ratio::Float64 #τ
  transfo_phase::Float64 #θ
  power_magnitude_max::Float64 #Smax
end


"""
    create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String) where T <: MatpowerLine_π_Smax

Create voltage variables at origin and destination of `link` in `link_vars`. \n
Return nothing.
"""
function create_vars!(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String) where T <: MatpowerLine_π_Smax
  origid, destid = link.orig, link.dest
  link_vars["Volt_orig"] = Variable(variable_name("VOLT", origid, "", scenario), Complex)
  link_vars["Volt_dest"] = Variable(variable_name("VOLT", destid, "", scenario), Complex)
  return
end


"""
    Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T <: MatpowerLine_π_Smax

Return the polynomial power at the origin of the line `link` depending quadratically on the voltages :
```julia
τ = element.transfo_ratio
θ = element.transfo_phase *pi/180
rs = element.resistance
xs = element.reactance
ys = 1/(rs+im*xs)
Yff = (ys+im*bc/2)/(τ^2)
Yft = -ys*1/(τ*exp(-im*θ))
Sorig = Yff*link_variables["Volt_orig"] + Yft*link_variables["Volt_dest"]
```
"""
function Sorig(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T <: MatpowerLine_π_Smax
  return element.baseMVA * link_vars["Volt_orig"] * conj(get_IorigMP(element, link_vars))
end

"""
    Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T <: MatpowerLine_π_Smax

Return the polynomial power at the destination of the line `link` depending quadratically on the voltages :
```julia
τ = element.transfo_ratio
θ = element.transfo_phase *pi/180
rs = element.resistance
xs = element.reactance
ys = 1/(rs+im*xs)
Ytf = -ys*1/(τ*exp(im*θ))
Ytt = ys+im*bc/2
Sdest = Ytf*link_variables["Volt_orig"] + Ytt*link_variables["Volt_dest"]
```
"""
function Sdest(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}) where T <: MatpowerLine_π_Smax
  return element.baseMVA * link_vars["Volt_dest"] * conj(get_IdestMP(element, link_vars))
end


## 3. Constraints creation
# #power constraints
function constraint(element::T, link::Link, elemid::String, elem_formulation::Symbol, link_vars::SortedDict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: MatpowerLine_π_Smax
  cstrs = SortedDict{String, Constraint}()

  if element.power_magnitude_max != 0
    Sor = Sorig(element, link, elemid, elem_formulation, link_vars)
    Sde = Sdest(element, link, elemid, elem_formulation, link_vars)
    Smax = element.power_magnitude_max

    cstrs[get_Smax_orig_cstrname()] = (abs2(real(Sor)) + abs2(imag(Sor))) << Smax^2
    cstrs[get_Smax_dest_cstrname()] = (abs2(real(Sde)) + abs2(imag(Sde))) << Smax^2
  end  

  return cstrs
end



## Utils functions
function get_IorigMP(element::MatpowerLine_π_Smax, link_variables::SortedDict{String, Variable})
  τ = element.transfo_ratio
  θ = element.transfo_phase *pi/180
  if τ == 0
    τ = 1
  end

  rs = element.resistance
  xs = element.reactance
  bc = element.susceptance
  ys = 1/(rs+im*xs)
  Yff = (ys+im*bc/2)/(τ^2)
  Yft = -ys*1/(τ*exp(-im*θ))
  return Yff*link_variables["Volt_orig"] + Yft*link_variables["Volt_dest"]
end

function get_IdestMP(element::MatpowerLine_π_Smax, link_variables::SortedDict{String, Variable})
  τ = element.transfo_ratio
  θ = element.transfo_phase *pi/180
  if τ == 0
    τ = 1
  end

  rs = element.resistance
  xs = element.reactance
  bc = element.susceptance
  ys = 1/(rs+im*xs)
  Ytf = -ys*1/(τ*exp(im*θ))
  Ytt = ys+im*bc/2
  return Ytf*link_variables["Volt_orig"] + Ytt*link_variables["Volt_dest"]
end
