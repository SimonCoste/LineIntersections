
function draw(l::Line)
    x = collect(range(-1, stop=1, length=100))
    plot!(x, point_in_line(l, x))
end

function draw_lines(H...)
    r = 1.1
    c = palette(:Set1_8)

    x = collect(range(-r, stop=r, length=2))
    v = collect(range(0, stop=2*pi, length=100))
    dessin = plot(sin.(v), cos.(v), color=:gray) #plot the circle
    for (i, arg) in enumerate(H)
        for l in arg
            plot!(dessin, x, point_in_line(l, x), color=c[i%8])
        end
    end
    plot!(dessin, xlims=(-r, r), ylims=(-r, r), legend=:none, aspect_ratio=:equal)
    return dessin
end

