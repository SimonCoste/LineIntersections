using DataStructures, Plots
#include("generator.jl")
#include("line_tools.jl")

import Base:<

global 𝐱 = -2.0

function <(
    ℓ₁::Line, 
    ℓ₂::Line
    )
    """
    Dynamic order on Lines, depending on the float 𝐱. 
    ℓ₁ < ℓ₂ is true when the intersection of ℓ₁ with
    the vertical line X with abscisse 𝐱 is below the intersection of ℓ₂
    with X. 
    """
    α = point_in_line(ℓ₁, 𝐱)
    β = point_in_line(ℓ₂, 𝐱)
    return α<β
end


function check_intersection_and_insert!(X, 
    X_dict, 
    intersections, 
    l1::Line, 
    l2::Line,
    r)
    """
    Tests if two lines intersect in the ball, then adds them to X and X_dict.
    """
    point = l1 ∩ l2
    if normof(point)<r
        insert!(X, point)
        X_dict[point] = Set([l1, l2])
        push!(intersections, point)
    end
end

function handle_intersection_point!(L, 
    X, 
    X_dict, 
    intersections, 
    ℓ1, 
    ℓ2, 
    absc)
    """
    This function handles the case of intersection points. 
    In particular it switches the lines who intersect by setting the sweep line 
    slightly to the right of its real position, but to the left of the next event
    which will occur.  
    """
    previous_rank = [sorted_rank(L, ℓ1), sorted_rank(L, ℓ2)]

    if previous_rank[1]<previous_rank[2]
        below = ℓ2
        above = ℓ1
    else
        below = ℓ1
        above = ℓ2
    end

    up_index = maximum(previous_rank)
    bel_index = minimum(previous_rank)
    
    if bel_index!=1
        line_to_test_below = L[bel_index-1]
        intersection_point = line_to_test_below ∩ below
        if normof(intersection_point)<1 && intersection_point[1] > absc
            if !haskey(X, intersection_point)
                insert!(X, intersection_point)
                X_dict[intersection_point]=Set([below, line_to_test_below])
                push!(intersections, intersection_point)
            end
        end
    end

    if up_index!=L.count 
        line_to_test_above = L[up_index + 1]
        intersection_point = line_to_test_above ∩ above
        if normof(intersection_point)<1 && intersection_point[1] > absc
            if !haskey(X, intersection_point)
                insert!(X, intersection_point)
                X_dict[intersection_point]=Set([above, line_to_test_above])
                push!(intersections, intersection_point)
            end
        end
    end

    delete!(L, ℓ1)
    delete!(L, ℓ2)
    
    next_future_event = X[1][1] 
    global 𝐱 = 0.5 * (absc + next_future_event)
    
    insert!(L, ℓ1)
    insert!(L, ℓ2)
end

sort_by_first_coordinate(a::Array{Float64, 2}) = a[:, sortperm(a[1, :])]


function intersections_bentley_ottman(H::Array{Line} ; r=1.)
    """
    A first (naive) implementation of the bentley_ottman algorithm for 
    computing the intersections in a circle of radius r. 

    Inputs:
        - H = an array or a list of Lines
        - r, a positive number (default 1).

    Outputs:
        - a 2×M array whose columns are the coordinates of the M intersection points of the lines in H.

    Current version does not support multiple intersection points nor vertical lines. 
    """


    ℒ = AVLTree{Line}() #Y-structure
    X = AVLTree{Array{Float64}}()
    X_dict = Dict{Array{Float64, 1}, Union{Int, Set{Line}}}()
    intersections  = Array{Array{Float64, 1}, 1}()
    step = 0
    
    points = get_endpoints(H)
    for i in 1:2*size(H)[1]
        p = points[:, i]
        insert!(X, p)
        X_dict[p] = i
    end


    while X.count > 0
        step = step + 1
        p = X[1]
        delete!(X, p)
        i = X_dict[p]

        if typeof(i)==Set{Line} #intersection point
            
            L = Array{Line}(undef, 2)
            k=1
            for line in i
                L[k] = line
                k = k+1
            end

            handle_intersection_point!(ℒ, X, X_dict, intersections, L[1], L[2], p[1])
        elseif i%2!=1 #endpoint
            
            ℓ = H[(i+1)÷2]
            line_index = sorted_rank(ℒ, ℓ)
            global 𝐱 = p[1]

            if line_index!=1 && line_index!=ℒ.count
                successor = ℒ[line_index+1]
                predecessor = ℒ[line_index-1]
                delete!(ℒ, ℓ)
                intersection_point = successor ∩ predecessor
                if normof(intersection_point)<1 && intersection_point[1]>𝐱
                    if !haskey(X, intersection_point)
                        insert!(X, intersection_point)
                        X_dict[intersection_point]=Set([successor, predecessor])
                        push!(intersections, intersection_point)
                    end
                end

            else
                delete!(ℒ, ℓ)
            end

        else #beginning point

            global 𝐱 = p[1]

            ℓ = H[(i+1)÷2]
            insert!(ℒ, ℓ)
            line_index = sorted_rank(ℒ, ℓ)
            
            if ℒ.count>1
                if line_index==1 
                    successor = ℒ[line_index+1]
                    check_intersection_and_insert!(X, X_dict, intersections, ℓ, successor, r)
                elseif line_index==ℒ.count
                    predecessor = ℒ[line_index-1]
                    check_intersection_and_insert!(X, X_dict, intersections, ℓ, predecessor, r)
                else
                    successor = ℒ[line_index+1]
                    predecessor = ℒ[line_index-1]
                    check_intersection_and_insert!(X, X_dict, intersections, ℓ, successor, r)
                    check_intersection_and_insert!(X, X_dict, intersections, ℓ, predecessor, r)
                end
            end
        end

    end
    
    res = hcat(intersections...)
    return sort_by_first_coordinate(res)

end




function intersections_naive(H::Array{Line}; r=1.)
    s = size(H)[1]
    k = Int(s*(s-1)/2)
    res = Array{Array{Float64, 1}, 1}()
    for i in 1:s
        for j in 1:(i-1)
            p = H[i] ∩ H[j]
            if normof(p)<r
                push!(res, p)
            end
        end
    end
    res = hcat(res...)
    return sort_by_first_coordinate(res)
end







#n = 10 #number of lines
#H = [Line(pi/4, 0.1), Line(-pi/2, 0.9)]
#H = hyperplanes_poisson(1, n)
#H = [ Line(4.0461103038672945, 0.7829732875581323), Line(1.8885159415469313, 0.455802052465595), Line(3.4109485517087896, 0.16889664754111178)]
#H = hyperplanes_poisson(1, n)


#bores = intersections_bentley_ottman(H)
#naiveres = intersections_naive(H)

# dessin = draw_lines(H)
# points = hcat(collect(keys(res))...)
# scatter!(points[1,:], points[2,:], 
# markersize=0.9)
# plot!(legend=:none, xlims=(-1, 1), ylims=(-1,1),
# aspect_ratio=:equal)
# dessin