"""
    mutable struct RTEMatpowerGenerator <: AbstractNodeLabel

Mutable structure descendant of AbstractNodeLabel

# Fields
- `power_min::Complex128`
- `power_max::Complex128`
- `cost_degree::Int64`
- `cost_coeffs::Array{Float64}`
- `gen_status::Bool`: generator on/off

"""
mutable struct RTEMatpowerGenerator <: AbstractNodeLabel
  power_min::Complex #Smin = Pmin + i Qmin
  power_max::Complex #Smax = Pmax + i Qmax
  active_power_input::Float64
  gen_status::Bool
end


"""
    create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String) where T <: RTEMatpowerGenerator

Create power variables n for generator `elemid` at `bus` in `bus_vars` if `elem_formulation == :NewVar`\n
Return nothing

"""
function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String) where T <: RTEMatpowerGenerator
    bus_vars[elemid] = Variable(variable_name("Sgen", bus, elemid, scenario), Complex)
    #bus_vars["lambda_$(sign)"] = Variable("lambda_$(sign)", Real)
  return
end

"""
[cstrname, polynom, lb, ub] = Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T <: RTEMatpowerGenerator

Return the polynomial contribution in power of generator `elemid` at `bus`(name, value, lower bound, upper bound). Will be used to construct power balance constraints in polynomial problem.\n
If `elem_formulation == :NbMinVar`, return generator bounds `["UNIT", Polynomial() ,Smin, Smax]`\n
If `elem_formulation == :NewVar`, return -Sgen `["UNIT", -Sgen, 0, 0]`\n
Return no contribution `["", Polynomial(), 0, 0]` otherwise
"""
function Snodal(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}) where T <: RTEMatpowerGenerator
  cstrname = "UNIT"
  return [cstrname, -bus_vars[elemid], 0, 0]

end


"""
    constraint(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: RTEMatpowerGenerator

Return all the constraints defined by generator `elemid` at `bus`. Will be used to construct constraints in polynomial problem.\n
If `elem_formulation == :NewVar`, return generator bounds : "Genbounds" => Smin <= Sgen <= Smax\n
Return empty dictionary otherwise
"""
function constraint(element::T, busid::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String, OPFpbs::OPFProblems) where T <: RTEMatpowerGenerator
  #bounds on the reactive power
    Qmin = imag(element.power_min)
    Qmax = imag(element.power_max)
    reactive_power = 0.5*(-im)*(bus_vars[elemid] - conj(bus_vars[elemid]))
    Pmin = real(element.power_min)
    Pmax = real(element.power_max)
    Pinput = element.active_power_input
    if !element.gen_status
      Pmin = Pmax = 0
      Qmin = Qmax = 0
    end
    active_power = 0.5*(bus_vars[elemid]+conj(bus_vars[elemid]))
    return SortedDict{String, Constraint}( "Reactivepower" => Qmin << reactive_power << Qmax)

end


"""
    cost(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_elems_vars::SortedDict{String, Variable}, Snode::Polynomial, lb, ub) where T <: RTEMatpowerGenerator

Return the polynomial contribution in objective generator `elemid` at `bus`.
"""
function cost(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_elems_var::SortedDict{String, Variable}, Snode::Polynomial, lb, ub) where T <: RTEMatpowerGenerator
  Sgen = bus_elems_var[elemid]
  Pgen = 0.5*(Sgen+conj(Sgen))
  if !element.gen_status
    Pgen = 0
  end
  return Pgen
end
