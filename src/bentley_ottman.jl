using DataStructures, Plots
#include("generator.jl")
#include("line_tools.jl")

#########################################################################
# ------------------- TYPES FOR THE BO ALGORITHM  -----------------------
#########################################################################

import Base:<

struct Event
    point::Array{Float64, 1}
    nature::String
    assets::Union{RootedLine, Tuple{RootedLine, RootedLine}}
end

function <(e::Event, f::Event)
    return e.point < f.point
end

#########################################################################
# ------------------- PREPROCESSING  ------------------------------------
#########################################################################

function insert_endpoints_of_one_line!(X::AVLTree{Event}, l::RootedLine, r::Real)
    (startpoint, endpoint) = intersection_with_circle(l, r)
    startevent = Event(startpoint, "startpoint", l)
    endevent = Event(endpoint, "endpoint", l)
    insert!(X, startevent) ; insert!(X, endevent)
end

function initialize_event_queue(H::Array{Line, 1}, r::Real)
    X = AVLTree{Event}()
    t = SweepLine(-r-1)
    R = [RootedLine(l, t) for l in H]
    for l in R
        insert_endpoints_of_one_line!(X, l, r)
    end
    return R, X, t
end

#########################################################################
# ------------------- HANDLING EVENTS  ----------------------------------
#########################################################################

function check_intersection_and_update!(a::Real, X::AVLTree{Event}, list, l1::RootedLine, l2::RootedLine, r::Real)
    """
    Checks if l1 and l2 intersect to the right of the sweep line at t and inside the disk of 
    radius r. If yes, checks if the intersection point had already been seen. If no, 
    registers this point in the list and add the event to the X-structure. 
    """
    if !haskey(list, Set([l1, l2])) #do nothing
        point::Array{Float64} = l1 ∩ l2
        if normof(point) < r && point[1] > a
            parents = point_in_line(l1, a)<point_in_line(l2, a) ? (l1, l2) : (l2, l1)
            e = Event(point, "intersection", parents) 
            insert!(X, e) ; list[Set([l1, l2])] = point 
        end
    end
end

function handle_startpoint!(t::SweepLine, X::AVLTree{Event}, ℒ::AVLTree{RootedLine}, list, e::Event, r::Real)
    point = e.point
    ℓ = e.assets
    t.r = point[1] #update the sweep line
    insert!(ℒ, ℓ)
    line_index = sorted_rank(ℒ, ℓ)
    
    if ℒ.count>1
        if line_index==1 
            successor = ℒ[line_index+1]
            check_intersection_and_update!(t.r, X, list, ℓ, successor, r)
        elseif line_index==ℒ.count
            predecessor = ℒ[line_index-1]
            check_intersection_and_update!(t.r, X, list, ℓ, predecessor, r)
        else
            successor = ℒ[line_index+1]
            predecessor = ℒ[line_index-1]
            check_intersection_and_update!(t.r, X, list, ℓ, successor, r)
            check_intersection_and_update!(t.r, X, list, ℓ, predecessor, r)
        end
    end
end

function handle_endpoint!(t::SweepLine, X::AVLTree{Event}, ℒ::AVLTree{RootedLine}, list, e::Event, r::Real)
    point = e.point
    ℓ = e.assets
    t.r = point[1]    #update the sweep line
    line_index = sorted_rank(ℒ, ℓ)

    if line_index!=1 && line_index!=ℒ.count
        successor = ℒ[line_index+1]
        predecessor = ℒ[line_index-1]
        check_intersection_and_update!(t.r, X, list, successor, predecessor, r)
    end
    delete!(ℒ, ℓ)
end

function handle_intersection!(t::SweepLine, X::AVLTree{Event}, ℒ::AVLTree{RootedLine}, list, e::Event, r::Real)
    point = e.point 
    (lower, upper) = e.assets
    (lower_rank, upper_rank) = [sorted_rank(ℒ, lower), sorted_rank(ℒ, upper)] #simplify ?

    if lower_rank!=1
        check_intersection_and_update!(point[1], X, list, upper, ℒ[lower_rank-1], r)
    end

    if upper_rank!=ℒ.count
        check_intersection_and_update!(point[1], X, list, lower, ℒ[upper_rank+1], r)
    end

    delete!(ℒ, lower) ; delete!(ℒ, upper)
    next_future_event = X[1]
    mid = 0.5 * (point[1] + next_future_event.point[1])
    t.r = mid
    insert!(ℒ, upper) ; insert!(ℒ, lower)

end

#########################################################################
# -------------------MAIN FUNCTION --------------------------------------
#########################################################################

function find_all_intersections(H::Array{Line} ; r::Real = 1.)

    @assert r>0 throw(ArgumentError("Radius r must me >0."))

    R, X, t = initialize_event_queue(H, r)
    list = Dict{Set{RootedLine}, Array{Float64}}()
    Y = AVLTree{RootedLine}()

    while X.count > 0

        event = X[1] 
        delete!(X, event)

        if event.nature == "startpoint"
            handle_startpoint!(t, X, Y, list, event, r)
        elseif event.nature == "endpoint"
            handle_endpoint!(t, X, Y, list, event, r)
        else 
            handle_intersection!(t, X, Y, list, event, r)
        end
    end

    return list
end

function intersections_naive(H::Array{Line}; r=1.)
    """
    Loops all over the pairs of distinct lines, 
    yielding a O(n^2) algorithm. 
    """
    s = size(H)[1]
    k = Int(s*(s-1)/2)
    res = Set{Array{Float64, 1}}()
    for i in 1:s
        for j in 1:(i-1)
            p = H[i] ∩ H[j]
            if normof(p)<r
                push!(res, p)
            end
        end
    end

    return res
end





# H = [Line(0.7544959394550933, 0.15515333752706018)
# Line(1.6832151948996787, 0.04401847448027585)
# Line(5.698726857470466, 0.07601306086085136)
# Line(4.232076668938423, 0.3089918936601266)]



# #find_all_intersections(H[1:2])

# H = hyperplanes_poisson(1, 5)

# include("generator.jl")
# @time a = find_all_intersections(H)
# @time b = intersections_naive(H)
# length(a)==length(b)