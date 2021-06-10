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
  flag::String
  voltage_input::Complex
end


"""
    create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String) where T <: RTEMatpowerGenerator

Create power variables n for generator `elemid` at `bus` in `bus_vars` if `elem_formulation == :NewVar`\n
Return nothing

"""
function create_vars!(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_vars::SortedDict{String, Variable}, scenario::String) where T <: RTEMatpowerGenerator
    bus_vars[elemid] = Variable(variable_name("Sgen", bus, elemid, scenario), Complex)
    sign = element.flag
    bus_vars["lambda_$(sign)"] = Variable("lambda_$(sign)", Real)
    bus_vars["Volt"] = Variable(variable_name("VOLT", bus, "", scenario), Complex)
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
    sign = element.flag
    λ = bus_vars["lambda_$(sign)"]
    if sign == "minus"
      return SortedDict{String, Constraint}( "Reactivepower" => Qmin << reactive_power << Qmax,
      "Activepower" => -active_power + Pmin+2*(Pinput-Pmin)*λ == 0,
      "ctr_lambda_$(sign)" => 0 << λ << 0.5)
    elseif sign == "plus"
      return SortedDict{String, Constraint}( "Reactivepower" => Qmin << reactive_power << Qmax,
      "Activepower" => -active_power + 2*Pinput-Pmax+2*(Pmax-Pinput)*λ == 0,
      "ctr_lambda_$(sign)"=> 0.5 << λ << 1)
    else
      exit("Problem with flag")
    end
end


"""
    cost(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_elems_vars::SortedDict{String, Variable}, Snode::Polynomial, lb, ub) where T <: RTEMatpowerGenerator

Return the polynomial contribution in objective generator `elemid` at `bus`.
"""
function cost(element::T, bus::String, elemid::String, elem_formulation::Symbol, bus_elems_var::SortedDict{String, Variable}, Snode::Polynomial, lb, ub) where T <: RTEMatpowerGenerator
  #autre option : minimisation ecart tensions
  v = bus_elems_var["Volt"]
  v_input = element.voltage_input
  return abs2(v-v_input)
end
